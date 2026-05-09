#!/usr/bin/env bash
# ============================================================================
# 04_tree_building/tree_building.sh
# ----------------------------------------------------------------------------
# Purpose:    Build species trees under four concatenation modelling strategies (LG-partitioned, GHOST, PMSF [LG+C20+FO+R], Dayhoff6) plus ASTRAL-IV / wASTRAL; compute gCF/sCF, AU/SH/KH topology test, RF/GRF distances, and per-tree summary statistics (DVMC, treeness, ABS)
# Inputs:     13-loci_concat/matrix{50,60,70,80,90}/faa/FcC_supermatrix.phy + partition; gene trees from step 03
# Outputs:    iqtree/{1_partition,2_ghost,4_pmsf,5_dayhoff6,6_astral}/matrix_<m>/<model>_<m>*.{treefile,cf.tree}; topology test report; rfdist matrix; per-tree CSV summaries
# Paper ref:  Methods 'Phylogenetic inference'; Fig 3, 4, 5; Tables 1, 2, 3; Suppl Figs 2-6
# Software:   IQ-TREE v3, ASTRAL-IV (ASTER), wASTRAL (ASTER), ASTRAL v5.7.8, gotree, readal/Phylogears, PhyKIT
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
ASTRAL_JAR="${ASTRAL_JAR:-/path/to/astral.5.7.8.jar}"
PHYLOGEARS_DIR="${PHYLOGEARS_DIR:-/path/to/Phylogears-master}"
cd "$PROJECT_ROOT"

#### BUSCO Tree Building

# Building concatenated supermatrix tree using IQ-TREE
mkdir -p iqtree

################
# 1. Partition using iqtree LG model
################
# Create directory
mkdir -p iqtree/1_partition

# Run it parallel for all matrix
parallel -j 3 '
  mkdir -p iqtree/1_partition/matrix_{};
  iqtree \
    -s 13-loci_concat/matrix{}/faa/FcC_supermatrix.phy \
    -p 13-loci_concat/matrix{}/faa/FcC_supermatrix_partition.txt \
    -m MFP \
    --mset LG \
    --msub nuclear \
    -rclusterf 10 \
    -B 1000 \
    --alrt 1000 \
    -T 32 \
    -pre iqtree/1_partition/matrix_{}/lg_{}
' ::: 50 60 70 80 90

####################
# 2. GHOST (General Heterogeneous evolution On a Single Topology)
####################
# Create directory
mkdir -p iqtree/2_ghost
# Run it parallel for all matrix
parallel -j 3 '
  mkdir -p iqtree/2_ghost/matrix_{};
  iqtree \
    -s 13-loci_concat/matrix{}/faa/FcC_supermatrix.phy \
    -m LG+FO+H4 \
    -B 1000 \
    --alrt 1000 \
    -T 32 \
    -pre iqtree/2_ghost/matrix_{}/ghost_{}
' ::: 50 60 70 80 90

###################
# 3. posterior mean site frequency (PMSF)
###################
# Create directory
mkdir -p iqtree/4_pmsf
# Using parallel
parallel -j 2 '
  mkdir -p iqtree/4_pmsf/matrix_{};
  iqtree \
    -s 13-loci_concat/matrix{}/faa/FcC_supermatrix.phy \
    -ft iqtree/1_partition/matrix_{}/lg_{}.treefile \
    -n 0 \
    -m LG+C20+F+R \
    -T 32 \
    -pre iqtree/4_pmsf/matrix_{}/{}
' ::: 50 60 70 80 90

#
parallel -j 2 '
  mkdir -p iqtree/4_pmsf/matrix_{};
  iqtree \
    -s 13-loci_concat/matrix{}/faa/FcC_supermatrix.phy \
    -fs iqtree/4_pmsf/matrix_{}/{}.sitefreq \
    -m LG+C20+FO+R \
    -T 32 \
    -B 1000 \
    --alrt 1000 \
    -pre iqtree/4_pmsf/matrix_{}/pmsf_{}
' ::: 50 60 70 80 90

###################
# 4. Dayhoff6-recoding
###################
# Create directory
mkdir -p iqtree/5_dayhoff6

# Using fas, convert .phy into .fas first
parallel -j 3 '
    mkdir -p iqtree/5_dayhoff6/matrix_{};
    readal \
      -in 13-loci_concat/matrix{}/faa/FcC_supermatrix.phy \
      -out iqtree/5_dayhoff6/matrix_{}/FcC_supermatrix.fas \
      -fasta \
    && \
    perl ${PHYLOGEARS_DIR}/pgrecodeseq.pl \
         --type=ANY "ARNDCQEGHILKMFPSTWYVX-01223220144145000554?" \
         iqtree/5_dayhoff6/matrix_{}/FcC_supermatrix.fas \
         iqtree/5_dayhoff6/matrix_{}/Dayhoff6_{}.fas \
    && \
    iqtree -s iqtree/5_dayhoff6/matrix_{}/Dayhoff6_{}.fas \
           -m GTR+R \
           -B 1000 \
           --alrt 1000 \
           -T 32 \
           -pre iqtree/5_dayhoff6/matrix_{}/dayhoff6_{}
' ::: 50 60 70 80 90

################
# 5. gene tree conflicts (gCF/sCF)
################
parallel -j 3 '
  set -euo pipefail
  m="{}"
  matdir="bad/13-loci_concat/matrix${m}"
  list="${matdir}/busco${m}.loci.list"
  out="${matdir}/all.gene.tre"

  # sanity checks
  [[ -s "$list" ]] || { echo "[M${m}] ERROR: missing $list" >&2; exit 1; }

  : > "$out"
  # concatenate the first matching gene tree for each locus
  grep -v "^[[:space:]]*$" "$list" | while IFS= read -r loci; do
    [[ -z "$loci" ]] && continue
    cat "9-gene_trees/iqtree/genetree/${loci}.treefile" >> "$out"
  done
  echo "[M${m}] Wrote $(grep -c \";\" \"$out\" 2>/dev/null || echo 0) trees -> $out"
' ::: 50 60 70 80 90

parallel -j 5 '
  dir=iqtree/6_astral/matrix_{};
  iqtree \
      -t  "$dir/astral_SU_{}.tre" \
      --gcf 13-loci_concat/matrix{}/all.gene.tre \
      -s  13-loci_concat/matrix{}/faa/original_FcC_supermatrix.phy \
      -p  13-loci_concat/matrix{}/faa/FcC_supermatrix_partition.txt \
      --scf 100 \
      -T 20 \
      --prefix "$dir/original_astral_SU_{}_gcf"
' ::: 50 60 70 80 90

##############
# 6. Topology AU/SH/KH test
##################

# ===============================================================
# CONFIGURATION
# We MUST use Matrix 60 because we have the matching .sitefreq file.
# ===============================================================
MAT="60"

# Paths
ALIGN="13-loci_concat/matrix${MAT}/faa/original_FcC_supermatrix.phy"

if [[ -f "iqtree/4_pmsf/matrix_${MAT}/pmsf_${MAT}.sitefreq" ]]; then
    SITEFREQ="iqtree/4_pmsf/matrix_${MAT}/pmsf_${MAT}.sitefreq"
elif [[ -f "iqtree/4_pmsf/matrix_${MAT}/${MAT}.sitefreq" ]]; then
    SITEFREQ="iqtree/4_pmsf/matrix_${MAT}/${MAT}.sitefreq"
else
    echo "ERROR: Could not find sitefreq file for Matrix ${MAT}"
    exit 1
fi

OUT_DIR="iqtree/test/topology_test/global_test"
mkdir -p "$OUT_DIR"

# Output files
RAW_TREES="${OUT_DIR}/all_candidates_raw.trees"
MAP_FILE="${OUT_DIR}/all_candidates.map"

# ===============================================================
# STEP 1: Combine all trees using robust extraction
# ===============================================================
echo "Collecting all trees from source folders..."

# Truncate outputs
: > "$RAW_TREES"
: > "$MAP_FILE"

# Enable nullglob so patterns matching nothing don't break the loop
shopt -s nullglob

i=0

# Loop through all specific tree patterns
# Note: The wildcard 'matrix_*/pmsf_*_gcf.cf.tree' captures Matrix 60 automatically. Raw: All unfiltered datasets; Bad: Excluded datasets
for f in \
  iqtree/1_partition/matrix_*/original_lg_*.treefile \
  iqtree/2_ghost/matrix_*/original_ghost_*_gcf.cf.tree \
  iqtree/4_pmsf/matrix_*/original_pmsf_*_gcf.cf.tree \
  iqtree/5_dayhoff6/matrix_*/original_dayhoff6_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/original_astral_old_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/original_astral_SU_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral4/original_astral_IV_*_gcf.cf.tree \
  raw-loci_concat/lg/matrix_*/original_raw_lg_*_gcf.cf.tree \
  raw-loci_concat/ghost/matrix_*/original_raw_ghost_*_gcf.cf.tree \
  raw-loci_concat/pmsf/matrix_*/original_raw_pmsf_*_gcf.cf.tree \
  raw-loci_concat/dayhoff6/matrix_*/original_raw_dayhoff6_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral/original_raw_astral_old_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral4/original_raw_astral_IV_*_gcf.cf.tree \
  bad/lg/matrix_*/original_bad_lg_*_gcf.cf.tree \
  bad/ghost/matrix_*/original_bad_ghost_*_gcf.cf.tree \
  bad/pmsf/matrix_*/original_bad_pmsf_*_gcf.cf.tree \
  bad/dayhoff6/matrix_*/original_bad_dayhoff6_*_gcf.cf.tree \
  bad/astral/matrix_*/astral/original_bad_astral_old_*_gcf.cf.tree \
  bad/astral/matrix_*/astral4/original_bad_astral_IV_*_gcf.cf.tree
do
  # 1) Extract Newick string
  newick=$(awk '/^[[:space:]]*\(.*;[[:space:]]*$/{print; exit}' "$f")

  if [[ -z $newick ]]; then
    echo "WARN: no Newick line found in $f" >&2
    continue
  fi

  i=$((i+1))
  echo "$newick" >> "$RAW_TREES"

  # 2) Extract Matrix Number
  mat=$(echo "$f" | grep -oE 'matrix_[0-9]+' | head -n1 | sed 's/matrix_//')
  [[ -z $mat ]] && mat="NA"

  # 3) Build Model Code
  base=$(basename "$f")
  case "$base" in
    original_lg_*)           code="lg" ;;
    original_ghost_*)        code="ghost" ;;
    original_pmsf_*)         code="pmsf" ;;
    original_dayhoff6_*)     code="dayhoff6" ;;
    original_raw_lg_*)       code="raw_lg" ;;
    original_raw_ghost_*)    code="raw_ghost" ;;
    original_raw_pmsf_*)     code="raw_pmsf" ;;
    original_raw_dayhoff6_*) code="raw_dayhoff6" ;;
    original_bad_lg_*)       code="bad_lg" ;;
    original_bad_ghost_*)    code="bad_ghost" ;;
    original_bad_pmsf_*)     code="bad_pmsf" ;;
    original_bad_dayhoff6_*) code="bad_dayhoff6" ;;


    # Logic for ASTRAL
    original_astral_SU_*)    code="wastral" ;;
    original_astral_old_*) code="astral" ;;
    original_astral_IV_*)    code="astral4" ;;
    original_raw_astral_old_*)    code="raw_astral" ;;
    original_raw_astral_IV_*)    code="raw_astral4" ;;
    original_bad_astral_old_*)    code="bad_astral" ;;
    original_bad_astral_IV_*)    code="bad_astral4" ;;

    *)              code="unknown" ;;
  esac

  name="${code}_${mat}"

  # 4) Write to map
  printf "Tree_%03d\t%s\t%s\n" "$i" "$name" "$f" >> "$MAP_FILE"
done

# ===============================================================
# STEP 2: Run the Global AU Test
# ===============================================================

iqtree -s 13-loci_concat/matrix60/faa/original_FcC_supermatrix.phy \
       -m LG+C20+FO+R+PMSF \
       -fs iqtree/4_pmsf/matrix_60/60.sitefreq \
       -z iqtree/test/topology_test/global_test/all_candidates_raw.trees \
       -n 0 -zb 10000 -zw -au \
       -T 32 \
       --prefix iqtree/test/topology_test/global_test/global_matrix60

####################
# 7. Robinson-Foulds distance
######################
# Create output directory
mkdir -p bad/test/rf_distance_test_no_wastral

# Define output files
out="bad/test/rf_distance_test_no_wastral/all_trees.tree"
map="bad/test/rf_distance_test_no_wastral/all_trees.map"

# Truncate (empty) output files before starting
: > "$out"
: > "$map"

# Turn on nullglob: if a wildcard matches nothing, it won't be treated as a literal string
shopt -s nullglob

i=0

# Loop through all specific tree patterns
for f in \
  iqtree/1_partition/matrix_*/lg_*_gcf.cf.tree \
  iqtree/2_ghost/matrix_*/ghost_*_gcf.cf.tree \
  iqtree/4_pmsf/matrix_*/pmsf_*_gcf.cf.tree \
  iqtree/5_dayhoff6/matrix_*/dayhoff6_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/unroot_astral_old_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral/unroot_astral_SU_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral4/unroot_astral_IV_*_gcf.cf.tree \
  raw-loci_concat/lg/matrix_*/raw_lg_*_gcf.cf.tree \
  raw-loci_concat/ghost/matrix_*/raw_ghost_*_gcf.cf.tree \
  raw-loci_concat/pmsf/matrix_*/raw_pmsf_*_gcf.cf.tree \
  raw-loci_concat/dayhoff6/matrix_*/raw_dayhoff6_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral/unroot_raw_astral_old_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral4/unroot_raw_astral_IV_*_gcf.cf.tree \
  bad/lg/matrix_*/bad_lg_*_gcf.cf.tree \
  bad/ghost/matrix_*/bad_ghost_*_gcf.cf.tree \
  bad/pmsf/matrix_*/bad_pmsf_*_gcf.cf.tree \
  bad/dayhoff6/matrix_*/bad_dayhoff6_*_gcf.cf.tree \
  bad/astral/matrix_*/astral/unroot_bad_astral_old_*_gcf.cf.tree \
  bad/astral/matrix_*/astral4/unroot_bad_astral_IV_*_gcf.cf.tree
do
  # 1) Robustly extract the Newick string
  # tr -d '\n' removes newlines
  newick=$(cat "$f" | tr -d '\n')

  # 2) Check if it looks like a tree
  # FIX: We escaped the semicolon (\;) so Bash doesn't think the command ends there.
  if [[ ! "$newick" =~ ^\(.*\;$ ]]; then
    printf 'WARN: Valid Newick tree not found in %s\n' "$f" >&2
    continue
  fi

  # Increment counter
  i=$((i+1))
  echo "$newick" >> "$out"

  # 3) Extract Matrix Number
  mat=$(echo "$f" | grep -oE 'matrix_[0-9]+' | head -n1 | sed 's/matrix_//')

  if [[ -z $mat ]]; then
      mat="NA"
  fi

  # 4) Build the Model Code
  base=$(basename "$f")

  case "$base" in
    lg_*)           code="lg" ;;
    ghost_*)        code="ghost" ;;
    pmsf_*)         code="pmsf" ;;
    dayhoff6_*)     code="dayhoff6" ;;
    raw_lg_*)       code="raw_lg" ;;
    raw_ghost_*)    code="raw_ghost" ;;
    raw_pmsf_*)     code="raw_pmsf" ;;
    raw_dayhoff6_*) code="raw_dayhoff6" ;;
    bad_lg_*)       code="bad_lg" ;;
    bad_ghost_*)    code="bad_ghost" ;;
    bad_pmsf_*)     code="bad_pmsf" ;;
    bad_dayhoff6_*) code="bad_dayhoff6" ;;


    # Logic for ASTRAL
    unroot_astral_SU_*)    code="wastral" ;;
    unroot_astral_old_*) code="astral" ;;
    unroot_astral_IV_*)    code="astral4" ;;
    unroot_raw_astral_old_*)    code="raw_astral" ;;
    unroot_raw_astral_IV_*)    code="raw_astral4" ;;
    unroot_bad_astral_old_*)    code="bad_astral" ;;
    unroot_bad_astral_IV_*)    code="bad_astral4" ;;

    *)              code="unknown" ;;
  esac

  # Construct the final name
  name="${code}_${mat}"

  # 5) Write the map line
  printf "Tree_%03d\t%s\t%s\n" "$i" "$name" "$f" >> "$map"
done

# Run the rf test
iqtree -t bad/test/rf_distance_test_no_wastral/all_trees.tree \
       -rf_all \
       -T 30 \
       -pre bad/test/rf_distance_test_no_wastral/rf_all

# Plot the result
python plot_rf_distance.py \
  --rfdist bad/test/rf_distance_test_no_wastral/rf_all.rfdist \
  --map    bad/test/rf_distance_test_no_wastral/all_trees.map \
  --out    bad/test/rf_distance_test_no_wastral/rf_mds_named.png \
  --coords bad/test/rf_distance_test_no_wastral/rf_mds_coords.tsv
# Compute GRF distances in R (TreeDist), then plot in Python
Rscript grf_dis_plot_clean.R \
  --tree   bad/test/rf_distance_test_no_astral/all_trees.tree \
  --coords bad/test/rf_distance_test_no_astral/grf_mds_coords.tsv

python plot_grf_distance.py \
  --coords_in bad/test/rf_distance_test_no_astral/grf_mds_coords.tsv \
  --map       bad/test/rf_distance_test_no_astral/all_trees.map \
  --out       bad/test/rf_distance_test_no_astral/grf_mds_FINAL.png

#########################
# 8. ASTRAL species-tree estimation
#########################

# Make the gene tree files
parallel -j 3 '
  set -euo pipefail
  m="{}"
  matdir="13-loci_concat/matrix${m}"
  list="${matdir}/busco${m}.loci.list"
  out="${matdir}/all.gene.tre"

  # sanity checks
  [[ -s "$list" ]] || { echo "[M${m}] ERROR: missing $list" >&2; exit 1; }

  : > "$out"
  # concatenate the first matching gene tree for each locus
  grep -v "^[[:space:]]*$" "$list" | while IFS= read -r loci; do
    [[ -z "$loci" ]] && continue
    cat "9-gene_trees/iqtree/genetree/${loci}.treefile" >> "$out"
  done
  echo "[M${m}] Wrote $(grep -c \";\" \"$out\" 2>/dev/null || echo 0) trees -> $out"
' ::: 50 60 70 80 90

### Run the Weighted ASTRAL, Use supports in the gene trees (e.g., IQ-TREE -abayes or bootstraps)
parallel -j 5 '
  set -euo pipefail
  m="{}"
  alltre="13-loci_concat/matrix${m}/all.gene.tre"
  outdir="iqtree/6_astral/matrix_${m}"
  mkdir -p "$outdir"

  [[ -s "$alltre" ]] || { echo "[M${m}] ERROR: missing $alltre" >&2; exit 1; }

  wastral -i "$alltre" -B -o "${outdir}/astral_${m}.tre" 2> "${outdir}/wastral_${m}.log"
  echo "[M${m}] wASTRAL done -> ${outdir}/astral_${m}.tre"
' ::: 50 60 70 80 90

# Add SU (substitutions/site) branch lengths to an existing topology
parallel -j 3 '
  set -euo pipefail
  m="{}"
  alltre="13-loci_concat/matrix${m}/all.gene.tre"
  topo="iqtree/6_astral/matrix_${m}/astral_${m}.tre"
  out="iqtree/6_astral/matrix_${m}/astral_SU_${m}.tre"
  log="iqtree/6_astral/matrix_${m}/castles_${m}.log"

  [[ -s "$topo"  &&  -s "$alltre" ]] || { echo "[M${m}] ERROR: need $topo and $alltre" >&2; exit 1; }

  # -C enables CASTLES-II; -c uses the existing topology
  astral4 -C -c "$topo" -i "$alltre" -o "$out" 2> "$log"
  echo "[M${m}] CASTLES-II done -> $out"
' ::: 50 60 70 80 90

### 1. Infer Topology using ASTRAL-IV (Unweighted)
# This replicates the robust logic of "Old ASTRAL" but is faster.
parallel -j 5 '
  set -euo pipefail
  m="{}"
  alltre="13-loci_concat/matrix${m}/all.gene.tre"
  outdir="iqtree/6_astral/matrix_${m}/astral4"
  mkdir -p "$outdir"

  # Run standard ASTRAL-IV
  # Note: No "-C" here. Just topology inference.
  astral4 -i "$alltre" -o "${outdir}/astral_${m}.tre" 2> "${outdir}/astral_${m}.log"

  echo "[M${m}] ASTRAL-IV Topology done -> ${outdir}/astral_${m}.tre"
' ::: 50 60 70 80 90

# 2. Calculate Branch Lengths (CASTLES-II) on that topology
# This step is correct in your previous script, just ensure it uses the input from step 1.
parallel -j 5 '
  set -euo pipefail
  m="{}"
  alltre="13-loci_concat/matrix${m}/all.gene.tre"
  # Use the topology we just made in step 1
  topo="iqtree/6_astral/matrix_${m}/astral4/astral_${m}.tre"
  out="iqtree/6_astral/matrix_${m}/astral4/astral_IV_${m}.tre"
  log="iqtree/6_astral/matrix_${m}/castles_${m}.log"

  # -C enables CASTLES, -c fixes the topology to the one provided
  astral4 -C -c "$topo" -i "$alltre" -o "$out" 2> "$log"
  echo "[M${m}] CASTLES-II Lengths done -> $out"
' ::: 50 60 70 80 90

# Run the old ASTRAL version
ASTRAL_JAR=${PROJECT_ROOT}/Astral/astral.5.7.8.jar

parallel -j 5 '
  set -euo pipefail
  alltre="13-loci_concat/matrix{}/all.gene.tre"
  outdir="iqtree/6_astral/matrix_{}"
  mkdir -p "$outdir"
  [[ -s "$alltre" ]] || { echo "[M{}] ERROR: missing $alltre" >&2; exit 1; }

  java -Xmx32g -jar "'"$ASTRAL_JAR"'" \
       -i "$alltre" \
       -o "${outdir}/astral_{}_old.tre" \
       2> "${outdir}/astral_{}.log"

  echo "[M{}] ASTRAL done -> ${outdir}/astral_{}_old.tre"
' ::: 50 60 70 80 90

# Unroot the ASTRAL trees

# Loop through all three tree file patterns
for file in iqtree/6_astral/matrix_*/astral_old_*_gcf.cf.tree \
            iqtree/6_astral/matrix_*/astral_SU_*_gcf.cf.tree \
            iqtree/6_astral/matrix_*/astral4/astral_IV_*_gcf.cf.tree; do

    # 1. Skip if file doesn't exist
    [ -e "$file" ] || continue

    # 2. Construct new filename (unroot_ at the start)
    dir=$(dirname "$file")      # Get the folder path
    base=$(basename "$file")    # Get the file name
    outfile="$dir/unroot_$base" # Combine them
    echo "  -> Output: $outfile"

    # 3. Run gotree unroot
    gotree unroot -i "$file" -o "$outfile"

done

########## 9. Run the ABS, GCF, SCF for all the trees
# Define Output File
OUTCSV="tree_support_stats.csv"

# Write Header
echo "TreeFile,Avg_SH_aLRT,Avg_UFBoot,Avg_gCF,Avg_sCF" > "$OUTCSV"

# Loop through all files matching your patterns
# We use 'ls' inside the loop to expand the wildcards (*)
for TREE in \
  iqtree/1_partition/matrix_*/lg_*_gcf.cf.tree \
  iqtree/2_ghost/matrix_*/ghost_*_gcf.cf.tree \
  iqtree/4_pmsf/matrix_*/pmsf_*_gcf.cf.tree \
  iqtree/5_dayhoff6/matrix_*/dayhoff6_*_gcf.cf.tree \
  raw-loci_concat/lg/matrix_*/raw_lg_*_gcf.cf.tree \
  raw-loci_concat/ghost/matrix_*/raw_ghost_*_gcf.cf.tree \
  raw-loci_concat/pmsf/matrix_*/raw_pmsf_*_gcf.cf.tree \
  raw-loci_concat/dayhoff6/matrix_*/raw_dayhoff6_*_gcf.cf.tree \
  bad/lg/matrix_*/bad_lg_*_gcf.cf.tree \
  bad/ghost/matrix_*/bad_ghost_*_gcf.cf.tree \
  bad/pmsf/matrix_*/bad_pmsf_*_gcf.cf.tree \
  bad/dayhoff6/matrix_*/bad_dayhoff6_*_gcf.cf.tree
do
  # Check if file exists (in case wildcard matches nothing)

  # 1. grep -o: extracts only strings looking like ")100/98/50.5/45.2"
  # 2. tr -d: removes the closing parenthesis ")"
  # 3. awk -F'/': splits by slash and calculates mean for columns 1,2,3,4

  VALS=$(grep -o ")[0-9.]\+/[0-9.]\+/[0-9.]\+/[0-9.]\+" "$TREE" \
    | tr -d ')' \
    | awk -F'/' '
      BEGIN { sum1=0; sum2=0; sum3=0; sum4=0; cnt=0 }
      {
        sum1+=$1;
        sum2+=$2;
        sum3+=$3;
        sum4+=$4;
        cnt++
      }
      END {
        if(cnt>0)
          printf "%.2f,%.2f,%.2f,%.2f", sum1/cnt, sum2/cnt, sum3/cnt, sum4/cnt;
        else
          print "NA,NA,NA,NA"
      }')

  echo "${TREE},${VALS}" >> "$OUTCSV"
done

# Run for astral trees
OUTCSV="astral_support_stats.csv"

# Write Header
# ASTRAL trees usually have: Posterior Probability (PP), gCF, sCF
echo "TreeFile,Avg_PP,Avg_gCF,Avg_sCF" > "$OUTCSV"

# Loop through all ASTRAL files
for TREE in \
  iqtree/6_astral/matrix_*/unroot_astral_old_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral4/unroot_astral_IV_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral_SU_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral/unroot_raw_astral_old_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral4/unroot_raw_astral_IV_*_gcf.cf.tree \
  bad/astral/matrix_*/astral/unroot_bad_astral_old_*_gcf.cf.tree \
  bad/astral/matrix_*/astral4/unroot_bad_astral_IV_*_gcf.cf.tree
do
  # Check if file exists

  # 1. grep -o: matches strings like ")1/32.2/41.7" (3 blocks of numbers)
  # 2. tr -d: removes the closing parenthesis ")"
  # 3. awk -F'/': splits by slash
  #    $1 = Posterior Probability (PP)
  #    $2 = gCF
  #    $3 = sCF

  VALS=$(grep -o ")[0-9.]\+/[0-9.]\+/[0-9.]\+" "$TREE" \
    | tr -d ')' \
    | awk -F'/' '
      BEGIN { sum_pp=0; sum_gcf=0; sum_scf=0; cnt=0 }
      {
        sum_pp+=$1;
        sum_gcf+=$2;
        sum_scf+=$3;
        cnt++
      }
      END {
        if(cnt>0)
          printf "%.4f,%.2f,%.2f", sum_pp/cnt, sum_gcf/cnt, sum_scf/cnt;
        else
          print "NA,NA,NA"
      }')

  echo "${TREE},${VALS}" >> "$OUTCSV"
done

############# 10. Run the DVMC for all the trees
# 1. Define Output
OUTCSV="tree_dvmc_stats.csv"
echo "TreeFile,DVMC" > "$OUTCSV"

# 2. Create a temporary file listing all the target trees
# We use 'ls' to expand the wildcards provided in your query
ls -1 \
  iqtree/1_partition/matrix_*/root_lg_*_gcf.cf.tree \
  iqtree/2_ghost/matrix_*/root_ghost_*_gcf.cf.tree \
  iqtree/4_pmsf/matrix_*/root_pmsf_*_gcf.cf.tree \
  iqtree/5_dayhoff6/matrix_*/root_dayhoff6_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/root_astral_old_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral_SU_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral4/root_astral_IV_*_gcf.cf.tree \
  raw-loci_concat/lg/matrix_*/root_raw_lg_*_gcf.cf.tree \
  raw-loci_concat/ghost/matrix_*/root_raw_ghost_*_gcf.cf.tree \
  raw-loci_concat/pmsf/matrix_*/root_raw_pmsf_*_gcf.cf.tree \
  raw-loci_concat/dayhoff6/matrix_*/root_raw_dayhoff6_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral/root_raw_astral_old_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral4/root_raw_astral_IV_*_gcf.cf.tree \
  bad/lg/matrix_*/root_bad_lg_*_gcf.cf.tree \
  bad/ghost/matrix_*/root_bad_ghost_*_gcf.cf.tree \
  bad/pmsf/matrix_*/root_bad_pmsf_*_gcf.cf.tree \
  bad/dayhoff6/matrix_*/root_bad_dayhoff6_*_gcf.cf.tree \
  bad/astral/matrix_*/astral/root_bad_astral_old_*_gcf.cf.tree \
  bad/astral/matrix_*/astral4/root_bad_astral_IV_*_gcf.cf.tree \
  2> /dev/null > tree_list.tmp

# 3. Run Parallel
# We use :::: to read from the temp file created above
parallel -j 20 '
  TREE={}

  # Calculate DVMC.
  # Note: DVMC requires the tree to be rooted to be meaningful.
  # If phykit fails or returns empty, we default to NA.

  DVMC=$(phykit degree_of_violation_of_a_molecular_clock "$TREE" 2>/dev/null)

  # Check if we got a number back
  if [[ -z "$DVMC" ]]; then
      DVMC="NA"
  fi

  echo "${TREE},${DVMC}"
' :::: tree_list.tmp >> "$OUTCSV"

# 4. Cleanup
rm tree_list.tmp

########## 11. Run the treeness for all the trees
# 1. Define Output
OUTCSV="tree_treeness_stats.csv"
echo "TreeFile,Treeness" > "$OUTCSV"

# 2. Create a temporary file listing all the target trees
# We use 'ls' to expand the wildcards provided in your query
ls -1 \
  iqtree/1_partition/matrix_*/lg_*_gcf.cf.tree \
  iqtree/2_ghost/matrix_*/ghost_*_gcf.cf.tree \
  iqtree/4_pmsf/matrix_*/pmsf_*_gcf.cf.tree \
  iqtree/5_dayhoff6/matrix_*/dayhoff6_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral_old_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/astral4/unroot_astral_IV_*_gcf.cf.tree \
  iqtree/6_astral/matrix_*/unroot_astral_SU_*_gcf.cf.tree \
  raw-loci_concat/lg/matrix_*/raw_lg_*_gcf.cf.tree \
  raw-loci_concat/ghost/matrix_*/raw_ghost_*_gcf.cf.tree \
  raw-loci_concat/pmsf/matrix_*/raw_pmsf_*_gcf.cf.tree \
  raw-loci_concat/dayhoff6/matrix_*/raw_dayhoff6_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral/raw_astral_old_*_gcf.cf.tree \
  raw-loci_concat/astral/matrix_*/astral4/unroot_raw_astral_IV_*_gcf.cf.tree \
  bad/lg/matrix_*/bad_lg_*_gcf.cf.tree \
  bad/ghost/matrix_*/bad_ghost_*_gcf.cf.tree \
  bad/pmsf/matrix_*/bad_pmsf_*_gcf.cf.tree \
  bad/dayhoff6/matrix_*/bad_dayhoff6_*_gcf.cf.tree \
  bad/astral/matrix_*/astral/bad_astral_old_*_gcf.cf.tree \
  bad/astral/matrix_*/astral4/unroot_bad_astral_IV_*_gcf.cf.tree \
  2> /dev/null > tree_list.tmp

# 3. Run Parallel
# Uses 'phykit treeness' to calculate the statistic
parallel -j 20 '
  TREE={}

  # Calculate Treeness
  # Returns a float (e.g., 0.65) representing sum of internal / sum of total branch lengths
  VAL=$(phykit treeness "$TREE" 2>/dev/null)

  # Check if we got a valid number back (handle potential errors)
  if [[ -z "$VAL" ]]; then
      VAL="NA"
  fi

  echo "${TREE},${VAL}"
' :::: tree_list.tmp >> "$OUTCSV"

# 4. Cleanup
rm tree_list.tmp

echo "Done. Results saved to $OUTCSV"