#!/usr/bin/env bash
# ============================================================================
# 03_extraction_filtering/BUSCO_extraction.sh
# ----------------------------------------------------------------------------
# Purpose:    Aggregate BUSCO loci, align (MAGUS), trim (trimAl), TreeShrink, gappy-seq filter, 8 sequence/tree quality filters, build per-locus gene trees, then concatenate into supermatrices at five taxon-occupancy thresholds (50-90%)
# Inputs:     busco_result/run_<species>/run_coleoptera_odb12/busco_sequences/*.faa; species.list
# Outputs:    13-loci_concat/matrix{50,60,70,80,90}/faa/FcC_supermatrix.phy + partition + .loci.list
# Paper ref:  Methods 'Locus filtering'; Fig 2; Tables 1, 3
# Software:   MAGUS, trimAl, IQ-TREE v3, TreeShrink, PhyKIT, csvtk, seqkit, FASconCAT-G
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
FASCONCAT_DIR="${FASCONCAT_DIR:-/path/to/FASconCAT-G}"
cd "$PROJECT_ROOT"

# BUSCO extraction and filtering (Single+Fragment BUSCO genes, no FNA)

MAFFT_DIR="$(dirname "$(which mafft)")"

################################
# 1) SPECIES & LOCI LISTS
################################
# Generate species list
ls busco_result/ > species.list
mapfile -t SPECIES_NAME < species.list

# Generate loci list using parallel processing
ROOT="busco_result"

awk -v R="$ROOT" 'NF{
  gsub(/\r$/,"")
  printf "%s/%s/run_coleoptera_odb12/busco_sequences\n", R, $0
}' species.list \
| xargs -I{} -P 8 bash -c '
  shopt -s nullglob
  for f in "$1"/*.faa; do
    bn=${f##*/}; printf "%s\n" "${bn%.*}"
  done
' _ {} \
| LC_ALL=C sort -u -o loci.list

mapfile -t LOCI_NAME < loci.list

################################
# 2) FASTA CONVERSION
################################
# Define a function for single-line FASTA + header renaming (parallel job)
#    - Convert multi-line FASTA to single-line
#    - Rename header to >SPECIES
#    - Remove '*' if the file is .faa (protein sequences)

convert_fasta() {
  file="$1"
  species="$(printf '%s\n' "$file" | sed -E 's|.*/(run_[^/]+)/run_coleoptera_odb12/.*|\1|')"

  awk -v S="$species" '
    BEGIN { RS=">"; ORS="" }
    NR>1 {
      n = split($0, a, /\r?\n/)        # n = number of lines in this record
      seq = ""
      for (i = 2; i <= n; i++) {       # join sequence lines
        gsub(/[ \t\r]/, "", a[i])
        seq = seq a[i]
      }
      gsub(/\*/, "", seq)              # drop stop codons
      if (length(seq) > 0) {           # **only print non-empty record**
        print ">" S "\n" seq "\n"
        exit                            # **first non-empty only**
      }
    }
  ' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
}

export -f convert_fasta

# Use parallel to process all .fna/.faa in single_copy_busco_sequences
#
#    -print0 / -0 safely handle filenames with spaces
#    -j CPU threads
# shows a progress bar

awk -v R="$ROOT" 'NF{
  gsub(/\r$/,""); s=$0; if (s !~ /^run_/) s="run_" s;
  printf "%s/%s/run_coleoptera_odb12/busco_sequences\n", R, s
}' species.list \
| parallel -j 16  --env convert_fasta '
    shopt -s nullglob
    for f in "{}"/*.faa; do
      convert_fasta "$f"
    done
'

################################
# 3) GATHER LOCI -> 1-raw_loci
################################
# Create output directories
mkdir -p 1-raw_loci/faa

# Create locus files in parallel
printf "%s\n" "${LOCI_NAME[@]}" | parallel -j 32 "touch 1-raw_loci/faa/{}.faa"

# Function to process each species
process_species() {
    SPECIES="$1"
    SC_DIR="busco_result/${SPECIES}/run_coleoptera_odb12/busco_sequences"
    touch 1-raw_loci/log.txt

    shopt -s nullglob   # patterns with no matches expand to nothing (prevents literal *.faa)
    for faa in "$SC_DIR"/*.faa; do
        locus="${faa##*/}"; locus="${locus%.faa}"
        cat "$faa" >> "1-raw_loci/faa/${locus}.faa"
    done
    shopt -u nullglob

    while read -r locus; do
        [[ ! -f "$SC_DIR/$locus.faa" ]] && echo "${locus} in ${SPECIES} does not exist" >> 1-raw_loci/log.txt
    done < loci.list
}

# Export the function for GNU Parallel
export -f process_species

# Run processing in parallel across species
parallel -j 30 process_species ::: "${SPECIES_NAME[@]}"

################################
# 4) FILTER LOCI
################################
#Filter loci having too few taxa (less than four)
mkdir -p 2-loci_filter/faa

# Count missing (by locus) from the log
awk 'NF{count[$1]++} END{for (l in count) print l, count[l]}' \
  1-raw_loci/log.txt > 2-loci_filter/missing_counts.tmp

# How many taxa (species) total?
TOTAL_TAXA=$(grep -v '^[[:space:]]*$' species.list | wc -l)
export TOTAL_TAXA

: > 2-loci_filter/sequence_number.log
: > 2-loci_filter/loci_name_filter.log

# Process each locus FASTA (faa only), in parallel
find 1-raw_loci/faa -type f -name '*.faa' -print0 \
| parallel -0 -j 20 '
    LOCI={/.}   # basename without .faa

    # Lookup missing count (default 0 if not found)
    MISSING=$(awk -v locus="$LOCI" '"'"'$1==locus {print $2}'"'"' 2-loci_filter/missing_counts.tmp)
    MISSING=${MISSING:-0}

    # Log the missing count
    printf "%s\t%s\n" "$LOCI" "$MISSING" >> 2-loci_filter/sequence_number.log

    # Keep loci present in ≥ (TOTAL_TAXA - 3) species
    if (( MISSING < TOTAL_TAXA - 3 )); then
      cp "1-raw_loci/faa/${LOCI}.faa" 2-loci_filter/faa/
      echo "$LOCI" >> 2-loci_filter/loci_name_filter.log
    fi
'

# Clean up
rm -f 2-loci_filter/missing_counts.tmp

################################
# 5) ALIGNMENT (16 cpu, 20 GB RAM per job)
################################
# Create output directory
mkdir -p 3-faa_align_magus
# Run MAGUS alignment in parallel
parallel -j 8 --colsep '\s+' '
  LOCUS={1}
  IN=2-loci_filter/faa/${LOCUS}.faa
  OUT=3-faa_align_magus/${LOCUS}.faa
  LOG=3-faa_align_magus/logs/${LOCUS}.log
  WORKDIR=3-faa_align_magus/workdir/${LOCUS}

  # 1) Clean any old workdir, then create it
  rm -rf "$WORKDIR"
  mkdir -p "$WORKDIR"

  # 2) Run MAGUS
  magus \
    -i "$IN" \
    -o "$OUT" \
    -np 16 \
    -d "$WORKDIR" \
    > "$LOG" 2>&1

  # 3) Check for success
  if [[ -s "$OUT" ]]; then
    rm -rf "$WORKDIR"
  else
  fi
' :::: 2-loci_filter/loci_name_filter.log

################################
# 6) TRIM
################################
##########Trim the alignments
mkdir -p 4-loci_trim/faa

# Read all filtered locus names into an array
mapfile -t LOCI_FILTER_ARR < 2-loci_filter/loci_name_filter.log

# Create a function for trimming one locus in parallel
trim_locus() {
    local locus="$1"
    trimal -in "3-faa_align_magus/${locus}.faa" \
                         -out "4-loci_trim/faa/${locus}.aa.fas" \
                         -automated1
    echo
}

export -f trim_locus

# Run trim jobs in parallel
parallel -j 50 trim_locus ::: "${LOCI_FILTER_ARR[@]}"

#########
# 7) Building gene trees (2 cpu, 1 GB RAM per job)
#############
# Create output directory for gene trees
# Run IQ-TREE for each locus in parallel
parallel -j 64 --colsep '\s+' '
  LOCUS={1}
  ALIGN=4-loci_trim/faa/${LOCUS}.aa.fas
  OUTDIR=5-gene_trees/iqtree
  PREFIX=${LOCUS}
  LOG=$OUTDIR/logs/${PREFIX}.log

  # Run IQ-TREE
  iqtree \
    -s $ALIGN \
    -m MFP \
    -mset Q.INSECT \
    -bb 1000 \
    -nt 2 \
    -pre $OUTDIR/${PREFIX} \
    > $LOG 2>&1 \
  && echo "Completed ${PREFIX}"
' :::: tree1.txt

# 1. Create the target directory (if it doesn’t already exist)
mkdir -p 5-gene_trees/iqtree/genetree

# 2. Move every .treefile into it
mv 5-gene_trees/iqtree/*.treefile 5-gene_trees/iqtree/genetree/

#################
# TREE SHRINK + Spurious homolog identification
##################
# Create output directory
mkdir -p 6-spurious_long_branch
#############
# 8) Run tree shrink for detecting abnormally long branches (0.5 GB per job)
############## TreeShrink for -m per-species ##################
OUTROOT=6-spurious_long_branch/1-treeshrink/per_species
LOCI=2-loci_filter/loci_name_filter.log
TREEDIR=5-gene_trees/iqtree/genetree
mkdir -p "$OUTROOT"

# ----------------------------------------
# 8.1) Build all.trees + index map (0-based)
#    Map format: <index>\t<locus>
# ----------------------------------------
MAP="${OUTROOT}/locus_index.tsv"
ALL="${OUTROOT}/all.trees"
: > "$MAP"
: > "$ALL"

awk -v d="$TREEDIR" -v map="$MAP" '
  BEGIN{ i=0 }
  NF>0 {
    loc=$0
    tf=sprintf("%s/%s.treefile", d, loc)
    # include only existing, non-empty trees (preserves order)
    cmd=sprintf("test -s \"%s\"", tf)
    if (system(cmd) == 0) {
      printf "%d\t%s\n", i, loc >> map   # record index -> locus
      print tf                           # <-- print paths to STDOUT (for the pipe)
      i++
    } else {
      printf "[WARN] missing or empty tree: %s\n", tf > "/dev/stderr"
    }
  }
' "$LOCI" \
| xargs -r cat > "$ALL"

TOTAL_TREES=$(wc -l < "$ALL" | tr -d ' ')
echo "[info] Wrote ${TOTAL_TREES} trees to $ALL"
echo "[info] Mapping saved to $MAP (index -> locus)"

# ----------------------------------------
# 8.2) Run TreeShrink per-species on ALL genes
#    Outputs: out_shrunk_<q>.trees + out_shrunk_RS_<q>.txt
# ----------------------------------------
THRESHOLDS=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

run_treeshrink.py \
  -t "$ALL" \
  -m per-species \
  -q "${THRESHOLDS[*]}" \
  -O out \
  -o "$OUTROOT" \
  > "${OUTROOT}/run.log" 2>&1

# ----------------------------------------
# 8.3) Per-threshold CSVs (Locus,SeqID) + summary
# ----------------------------------------
OUTROOT=6-spurious_long_branch/1-treeshrink/per_species
MAP="${OUTROOT}/locus_index.tsv"    # "<0-based-index>\t<locus>" from your build step
ALL="${OUTROOT}/all.trees"
[[ -s "$MAP" && -s "$ALL" ]] || { echo "Missing $MAP or $ALL"; exit 1; }

# Make a 1-based map for the line-per-tree format
MAP1="${OUTROOT}/locus_index_1based.tsv"
awk '{print NR "\t" $2}' "$MAP" > "$MAP1"

THRESHOLDS=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
TOTAL_LOCI=$(wc -l < "$ALL" | tr -d ' ')
echo -e "q\tremoved_seqs\taffected_loci\ttotal_loci" > "${OUTROOT}/summary.tsv"

# Helper: find the "removing set" file for a given q (new and old names)
rs_for_q () {
  local q="$1"
  local f
  for f in \
    "${OUTROOT}/out_${q}.txt" \
    "${OUTROOT}/out_$(printf '%.2f' "$q").txt" \
    "${OUTROOT}/out_$(printf '%.1f' "$q").txt" \
    "${OUTROOT}/out_shrunk_RS_${q}.txt" \
    "${OUTROOT}/out_shrunk_RS_$(printf '%.2f' "$q").txt" \
    "${OUTROOT}/out_shrunk_RS_$(printf '%.1f' "$q").txt"
  do
    [[ -s "$f" ]] && { echo "$f"; return 0; }
  done
  return 1
}

for q in "${THRESHOLDS[@]}"; do
  RS="$(rs_for_q "$q" || true)"
  OUTCSV="${OUTROOT}/treeshrink_${q}.csv"
  echo "Locus,SeqID" > "$OUTCSV"

  if [[ -z "${RS:-}" ]]; then
    printf "%s\t0\t0\t%s\n" "$q" "$TOTAL_LOCI" >> "${OUTROOT}/summary.tsv"
    echo "[warn] No removing-set file for q=${q}; header-only ${OUTCSV}"
    continue
  fi

  # Try FORMAT A: numeric table (older slides) → first token is an integer index
  HITS_A="${OUTROOT}/.hitsA_${q}.csv"
  awk 'NR==FNR {idx[$1]=$2; next}
       $1 ~ /^[0-9]+$/ {print idx[$1] "," $2}' "$MAP" "$RS" > "$HITS_A"

  # We parsed the table format successfully
  HITS="$HITS_A"

  # Loci with no removals at this q
  ALL_LOCI="${OUTROOT}/.all_${q}.txt"
  HIT_LOCI="${OUTROOT}/.hit_${q}.txt"
  NONE="${OUTROOT}/.none_${q}.csv"

  cut -f2 "$MAP" | sort > "$ALL_LOCI"
  cut -d, -f1 "$HITS" | sort -u > "$HIT_LOCI"
  comm -23 "$ALL_LOCI" "$HIT_LOCI" | awk '{print $0 ",None"}' > "$NONE"

  # Final CSV
  cat "$HITS" "$NONE" >> "$OUTCSV"

  removed_seqs=$(wc -l < "$HITS" | tr -d ' ')
  affected_loci=$(wc -l < "$HIT_LOCI" | tr -d ' ')
  printf "%s\t%s\t%s\t%s\n" "$q" "$removed_seqs" "$affected_loci" "$TOTAL_LOCI" >> "${OUTROOT}/summary.tsv"

  rm -f "$HITS_A" "$HITS_B" "$ALL_LOCI" "$HIT_LOCI" "$NONE"
  echo "[done] ${OUTCSV}"
done

#################
# 10) Filter out spurious homologs and long branches
####################
# Concatenate TreeShrink removals into the summary CSV
echo "Locus,SeqID" > 6-spurious_long_branch/summary_stats.csv

csvtk concat -U \
  6-spurious_long_branch/1-treeshrink/per_species/treeshrink_0.05.csv \
| csvtk uniq -H -f 1,2 \
>> 6-spurious_long_branch/summary_stats.csv

awk -F, 'BEGIN{OFS=","} NR==1 || $2!="None"' 6-spurious_long_branch/summary_stats.csv \
  > 6-spurious_long_branch/summary_stats_final.csv

# Create output dirs
mkdir -p \
  6-spurious_long_branch/faa

tail -n +2 6-spurious_long_branch/summary_stats_final.csv \
  | cut -d, -f1 \
  | sort -u \
  > 6-spurious_long_branch/loci_to_realign.txt

mkdir -p 6-spurious_long_branch/remove_lists

# 1) Read CSV and write SeqIDs to locus-specific files
awk -F, '
  NR>1 && $2!="None" {
    file = "6-spurious_long_branch/remove_lists/" $1 ".txt"
    print $2 >> file
    close(file)
  }
' 6-spurious_long_branch/summary_stats_final.csv

# Make sure top‐level dirs exist
mkdir -p tmp magus_work

# Parallel MAGUS realign
parallel -j 8 --colsep '\s+' '
  locus={1}
  rm_list=6-spurious_long_branch/remove_lists/${locus}.txt
  extract_fa=tmp/${locus}/${locus}.faa
  magus_out=tmp/${locus}/${locus}.magus.aa.fas
  workdir=magus_work/${locus}

  mkdir -p tmp/${locus} "$workdir"

  # remove unwanted seqs
  seqkit grep -n -v -f "$rm_list" 2-loci_filter/faa/${locus}.faa > "$extract_fa"

  # how many sequences left?
  n=$(grep -c "^>" "$extract_fa")

  if (( n < 3 )); then
      # not enough for MAGUS; just copy original/filtered file forward
      cp "$extract_fa" "$magus_out"
      echo "$locus – skipped MAGUS (only $n seqs)"
  else
      # realign with MAGUS
      magus -i "$extract_fa" -o "$magus_out" -np 16 -d "$workdir"
      echo "$locus – MAGUS realigned ($n seqs)"
  fi
' :::: 6-spurious_long_branch/loci_to_realign.txt

rm -rf magus_work

# Parallel trimAl
parallel -j 100 --colsep '\s+' '
  locus={}
  magus_aa=tmp/${locus}/${locus}.magus.aa.fas
  aa_out=6-spurious_long_branch/faa/${locus}.aa.fas

  [[ -e "$magus_aa" ]] || { echo "Missing input for $locus" >&2; exit 0; }

  # (1) trim columns once on the AA alignment
  trimal -in "$magus_aa" -out "$aa_out" -automated1

  rm -rf tmp/${locus} magus_work/${locus}
' :::: 6-spurious_long_branch/loci_to_realign.txt

rm -rf tmp

# Copy the rest of the loci
comm -23 <(sort -u 2-loci_filter/loci_name_filter.log) \
         <(sort -u 6-spurious_long_branch/loci_to_realign.txt) \
     > 6-spurious_long_branch/loci_to_copy.txt

parallel -j 100 '
  loc={}
  cp 4-loci_trim/faa/${loc}.aa.fas  6-spurious_long_branch/faa/
' :::: 6-spurious_long_branch/loci_to_copy.txt

##############################
# 11) Filter out the gappy alignments (70% coverage now)
##############################
mkdir -p 7-gappy_filt/faa

echo "Locus,SeqID" > 7-gappy_filt/removed_sequences_0.7.csv

parallel -j 100 --colsep '\s+' '
  LOCUS={1}
  PROT_IN=6-spurious_long_branch/faa/${LOCUS}.aa.fas
  PROT_OUT=7-gappy_filt/faa/${LOCUS}.aa.fas

  # run the filter script; append removed entries to the common CSV
  python filter_by_coverage.py \
    --in_prot "$PROT_IN" \
    --out_prot "$PROT_OUT" \
    --mincov 0.7 \
' :::: 2-loci_filter/loci_name_filter.log >> 7-gappy_filt/removed_sequences_0.7.csv

tail -n +2 7-gappy_filt/removed_sequences_0.7.csv \
  | cut -d, -f1 \
  | sort -u \
  > 7-gappy_filt/loci_to_filter.txt

# Parallel trimAl again
parallel -j 100 '
  locus={}

  # Current trimmed protein alignment
  aa_in=7-gappy_filt/faa/${locus}.aa.fas
  aa_tmp=7-gappy_filt/faa/${locus}.aa.trim.fas

  # Guard clause: need the AA file
  [[ -s "$aa_in" ]] || { echo "Missing AA input for ${locus}" >&2; exit 0; }

  # Trim the peptide alignment (heuristic method choice)
  trimal -in "$aa_in" -out "$aa_tmp" -automated1

  # Replace AA file with the newly trimmed one
  mv -f "$aa_tmp" "$aa_in"
' :::: 7-gappy_filt/loci_to_filter.txt

rm -rf tmp

##################
# 12) FILTER LOCI again
##################
#Filter loci having too few taxa (less than four)
mkdir -p 8-loci_filter/faa

# List number of sequences in each loci
find 7-gappy_filt/faa -name '*.aa.fas' -print0 |
  parallel -0 -j 100 -k '
    aa_file={}
    locus=$(basename "$aa_file" .aa.fas)

    # count headers (sequences)
    n_aa=$(grep -c "^>" "$aa_file")

    printf "%s\t%s\n" "$locus" "$n_aa"
  ' \
  > 8-loci_filter/sequence_number.txt

# Filter loci contains less than 4 sequences
find 7-gappy_filt/faa -name '*.aa.fas' -print0 |
  parallel -0 -j 100 '
    locus=$(basename {} .aa.fas)
    aa_in="7-gappy_filt/faa/${locus}.aa.fas"

    # count headers (sequences)
    n_aa=$(grep -c "^>" "$aa_in")

    if (( n_aa >= 4 )); then
        mkdir -p 8-loci_filter/faa
        cp "$aa_in" 8-loci_filter/faa/
        echo "Kept  $locus (AA:$n_aa)"
    else
        echo "Dropped $locus (AA:$n_aa)"
    fi
  '
# Write the filtered loci names to a log file
find 8-loci_filter/faa -name '*.aa.fas' -print0 \
| parallel -0 -j 100 'basename {} .aa.fas' \
| sort -u > 8-loci_filter/loci_name_filter_0.7.log

################
# 13) Building gene trees
################
# Create output directory for gene trees

cat 6-spurious_long_branch/loci_to_realign.txt \
    7-gappy_filt/loci_to_filter.txt       \
  | sort -u > 9-gene_trees/all_loci_to_fix.txt

comm -12 \
  <(sort -u 9-gene_trees/all_loci_to_fix.txt) \
  <(sort -u 8-loci_filter/loci_name_filter_0.7.log)  \
  > 9-gene_trees/loci_to_rerun.txt

# Run IQ-TREE for each locus in parallel (change to 2 cpu with 1 G RAM for final species)
parallel -j 64 --colsep '\s+' '
  LOCUS={1}
  ALIGN=8-loci_filter/faa/${LOCUS}.aa.fas
  OUTDIR=9-gene_trees/iqtree
  PREFIX=${LOCUS}
  LOG=$OUTDIR/logs/${PREFIX}.log

  # Run IQ-TREE
  iqtree \
    -s $ALIGN \
    -m MFP \
    -mset Q.INSECT \
    -bb 1000 \
    -nt 2 \
    -pre $OUTDIR/${PREFIX} \
    > $LOG 2>&1 \
  && echo "Completed ${PREFIX}"
' :::: 9-gene_trees/loci_to_rerun.txt

# 1. Create the target directory (if it doesn’t already exist)
mkdir -p 9-gene_trees/iqtree/genetree

# 2. Copy every .treefile into it
mv 9-gene_trees/iqtree/*.treefile 9-gene_trees/iqtree/genetree/

# Copy the previous treefile into it
comm -23 <(sort -u 8-loci_filter/loci_name_filter_0.7.log) \
         <(sort -u 9-gene_trees/loci_to_rerun.txt) \
  > 9-gene_trees/tree_to_copy.txt

parallel -j 100 '
  loc={}
  cp 5-gene_trees/iqtree/genetree/${loc}.treefile  9-gene_trees/iqtree/genetree/
' :::: 9-gene_trees/tree_to_copy.txt

###########
# Alignment based filtering
###############
# Create output directory for alignment filtering
mkdir -p 10-alnbased
##########
# 14)Filter by number of parsimony-informative sites (aa)
##########
# The number of parsimony informative sites in an alignment is associated with strong phylogenetic signal.
# Create output directory for parsimony filtering
mkdir -p 10-alnbased/1_parsimony_filter
# Calculate aligment length, num.pis and per.pis, integrate results to csv
OUTCSV=10-alnbased/1_parsimony_filter/parsimony_stats.csv

# Write header
echo "Locus,Parsimony_sites,Alignment_length,Percent_pis" > "$OUTCSV"
# Use phykit to calculate parsimony informative sites

parallel -j 120 '
  FASTA=8-loci_filter/faa/{}.aa.fas
  read PIS LEN PCT < <(phykit parsimony_informative_sites "$FASTA")
  echo "{},${PIS},${LEN},${PCT}"
' :::: 8-loci_filter/loci_name_filter_0.7.log \
  >> "$OUTCSV"

###############
# 15) Filter by Relative composition variability (RCV) (aa)
###############
# RCV describes the average variability in sequence composition among taxa. Lower RCV values are thought to be desirable because they represent a lower composition bias in an alignment.
# Create output directory for RCV filtering
mkdir -p 10-alnbased/2_RCV_filter
# Calculate RCV results to csv
OUTCSV=10-alnbased/2_RCV_filter/RCV_stats.csv

# Write header
echo "Locus,RCV" > "$OUTCSV"
# Use phykit

parallel -j 120 '
  L={}
  FASTA=8-loci_filter/faa/${L}.aa.fas

  # phykit rcv prints a single number: the relative composition variability
  RCV=$(phykit rcv "$FASTA")
  echo "${L},${RCV}"
' :::: 8-loci_filter/loci_name_filter_0.7.log \
  >> "$OUTCSV"

###############
# 16) Filter by symmetry tests against SRH hypotheses (aa)
###############
mkdir -p 10-alnbased/3_symtest
OUTCSV=10-alnbased/3_symtest/symtest_pvalues.csv

# Write header
echo "Locus,SRH_pval" > "$OUTCSV"

# Run per‐locus symmetry test in parallel
parallel -j 64 '
  LOCUS={1}
  ALN=8-loci_filter/faa/${LOCUS}.aa.fas
  PREFIX=10-alnbased/3_symtest/${LOCUS}
  LOG=${PREFIX}.log

  iqtree \
    -s "$ALN" \
    --symtest-only \
    -T 2 \
    -pre "$PREFIX" \
    >"$LOG" 2>&1 && echo "Completed ${LOCUS}"
' :::: 8-loci_filter/loci_name_filter_0.7.log

parallel -j 120 '
  f={}
  locus=$(basename "$f" .symtest.csv)
  # grab first real data line, cut out column 4
  pval=$(grep -vE "^(#|Name,)" "$f" | head -n1 | cut -d, -f4)
  printf "%s,%s\n" "$locus" "$pval"
' ::: 10-alnbased/3_symtest/*.symtest.csv \
>> "$OUTCSV"

##############
# 17) Filter by evolutionary rate (average pairwise identity) (aa)
##############
# Higher PI indicates slower evolution (more conserved sequences), whereas lower PI indicates faster evolution (more divergent sequences)
# Create output directory for evolutionary rate filtering
mkdir -p 10-alnbased/4_EvoRate_filter
# Calculate results to csv
OUTCSV=10-alnbased/4_EvoRate_filter/EvoRate_stats.csv

# Write header
echo "Locus,evolutionary_rate" > "$OUTCSV"
# Use phykit

parallel -j 120 '
  L={};
  FASTA=8-loci_filter/faa/${L}.aa.fas;
  MEAN=$(phykit pairwise_identity "$FASTA" \
           | awk "/^mean:/ { print \$2 }");
  printf "%s,%.4f\n" "$L" "$MEAN";
' :::: 8-loci_filter/loci_name_filter_0.7.log \
>> "$OUTCSV"

##############
# 18) Filter by Likelihood mapping (aa)
##############
# Higher value of Percent Fully Resolved Qurtet means stronger phylogenetic informativess/signal
# Create output directory
mkdir -p 10-alnbased/5_Likelihood

# Run Likelihood mapping in parallel (4-8G per job)
parallel -j 64 '
  LOCUS={1}
  ALN=8-loci_filter/faa/${LOCUS}.aa.fas
  PREFIX=10-alnbased/5_Likelihood/${LOCUS}
  LOG=${PREFIX}.log
  OUT_IQTREE=${PREFIX}.iqtree

  [[ -s "$ALN" ]] || { echo "[Warning] ${LOCUS} missing alignment" >&2; exit 0; }

  echo "Skip ${LOCUS} (exists: ${OUT_IQTREE})"
  exit 0

  mkdir -p "$(dirname "$PREFIX")"

  iqtree \
    -s "$ALN" \
    -mset Q.insect \
    -lmap 1000000 \
    -n 0 \
    -T 2 \
    -pre "$PREFIX" \
    -quiet \
    -redo \
    >"$LOG" 2>&1 && echo "Completed ${LOCUS}"
' :::: 8-loci_filter/loci_name_filter_0.7.log

comm -23 <(sort 8-loci_filter/loci_name_filter_0.7.log) \
         <(ls 10-alnbased/5_Likelihood/*.iqtree \
              | sed 's!.*/!!;s/\.iqtree$//' | sort) \
      > 10-alnbased/5_Likelihood/misstree.log

# Calculate results to csv
OUTCSV=10-alnbased/5_Likelihood/Likelihood_stats.csv
tmpdir=$(mktemp -d)

parallel -j 120 '
  fq={}; loc=${fq##*/}; loc=${loc%.iqtree}
  val="NA"
  if line=$(grep -m1 "Number of fully resolved  quartets" "$fq" 2>/dev/null); then
      [[ $line =~ ([0-9]+([.][0-9]+)?)\ *% ]] && val=${BASH_REMATCH[1]}
  fi
  printf "%s,%s\n" "$loc" "$val" > "'"$tmpdir"'/${loc}.csv"
' ::: 10-alnbased/5_Likelihood/*.iqtree

# assemble (header once)
{
  echo "Locus,Percent_FullyResolvedQuartet"
  cat "$tmpdir"/*.csv
} > "$OUTCSV"

rm -r "$tmpdir"

##############
# Tree-based filtering
##############
# Create output directory for tree-based filtering
mkdir -p 11-treebased

##################
# 19) Filter by average bootstraps (ABS)
###################
# average bootstraps (ABS) is the average of all bootstrap values in a tree. More values are thought to be of more phylogenetic signal.
mkdir -p 11-treebased/1_ABS_filter

OUTCSV=11-treebased/1_ABS_filter/ABS_stat.csv
echo "Locus,ABS" > "$OUTCSV"

parallel -j 120 '
  LOCUS={1}
  TREE="9-gene_trees/iqtree/genetree/${LOCUS}.treefile"

  # break on “)”, grab the bootstrap number before “:”, drop the first (root) line,
  # strip semicolons, then compute mean
  ABS=$(sed "s/)/\n/g" "$TREE" \
    | cut -d: -f1 \
    | sed "1d; s/;//g" \
    | awk '\''/^[0-9]/ { sum+=$1; cnt++ }
               END { if(cnt>0) printf("%.2f", sum/cnt); else print "NA" }'\'')

  echo "${LOCUS},${ABS}"
' :::: 8-loci_filter/loci_name_filter_0.7.log \
  >> "$OUTCSV"

##################
# 20) Degree of violation of the molecular clock (DVMC)
##################
# Degree of violation of the molecular clock (DVMC) describes degree of violation of the molecular clock in a phylogeny. Lower DVMC values are thought to be desirable because they represent a lower degree of violation of the molecular clock.
mkdir -p 11-treebased/2_DVMC_filter

OUTCSV=11-treebased/2_DVMC_filter/DVMC_stat.csv
echo "Locus,DVMC" > "$OUTCSV"

parallel -j 120 '
  L={}
  TREE="9-gene_trees/iqtree/genetree/${L}.treefile"

  # Compute DVMC using PhyKIT’s CLI
  DVMC=$(phykit degree_of_violation_of_a_molecular_clock "$TREE")

  # Emit one CSV row per locus
  echo "${L},${DVMC}"
' :::: 8-loci_filter/loci_name_filter_0.7.log \
  >> "$OUTCSV"

##################
# 21) Filter by Treeness: describes the proportion of the tree distance found on internal branches. Treeness can be used as a measure of the signal-to-noise ratio in a phylogeny.
##################
mkdir -p 11-treebased/3_treeness_filter

OUTCSV=11-treebased/3_treeness_filter/treeness_stat.csv
echo "Locus,treeness" > "$OUTCSV"

parallel -j 120 '
  L={}
  TREE="9-gene_trees/iqtree/genetree/${L}.treefile"

  # Compute treeness using PhyKIT’s CLI
  treeness=$(phykit treeness "$TREE")

  # Emit one CSV row per locus
  echo "${L},${treeness}"
' :::: 8-loci_filter/loci_name_filter_0.7.log \
  >> "$OUTCSV"

##################
# 22) Make a summary table of all align and tree-based filtering results
##################
# Create output directory for summary
mkdir -p 12-final-loci_filt
OUTCSV=12-final-loci_filt/summary_stats_tree.csv
touch "$OUTCSV"
# Combined all the stats into this summary table
csvtk join \
  -j 20 \
  -H \
  -f Locus \
  -o "$OUTCSV" \
  10-alnbased/1_parsimony_filter/parsimony_stats.csv \
  10-alnbased/2_RCV_filter/RCV_stats.csv \
  10-alnbased/3_symtest/symtest_pvalues.csv \
  10-alnbased/4_EvoRate_filter/EvoRate_stats.csv \
  10-alnbased/5_Likelihood/Likelihood_stats.csv \
  11-treebased/1_ABS_filter/ABS_stat.csv \
  11-treebased/2_DVMC_filter/DVMC_stat.csv \
  11-treebased/3_treeness_filter/treeness_stat.csv

awk -F, '
  NR==1 {print $0",treeness/RCV"; next}
  {print $0","$11/$5}   # replace 8 and 9 with the actual column numbers
' 12-final-loci_filt/summary_stats_tree.csv \
  > 12-final-loci_filt/summary_stats_tree_with_ratio.csv

##################
# 23) Filter loci
##################
# Create output directory for filtered loci
mkdir -p 12-final-loci_filt/faa
# Read the summary stats and filter based on criteria
# Filter criteria:
#   - Parsimony sites >= 100
#   - RCV: ≤0.2
#   - SRH p-value: >=0.05
#   - EvoRate: >0.5, <=0.9
#   - Likelihood: >= 50
#   - ABS >= 70
#   - DVMC: <= 0.4
#   - treeness/RCV: >= 0.7

# Make a list of loci that pass the criteria
OUTCSV=12-final-loci_filt/final.txt
touch "$OUTCSV"

# Filter based on the maintain.txt file
parallel -j 120 '
  cp "8-loci_filter/faa/{}.aa.fas" "12-final-loci_filt/faa/"
' :::: 12-final-loci_filt/final.txt

################################
# 24) FASconCAT: “ALL” DATASETS (FAA only) (20G RAM per job)
################################
export FASconCAT_DIR="${FASCONCAT_DIR}"

# 1. Setup directories
mkdir -p 13-loci_concat/all/faa
mkdir -p 13-loci_concat/matrix{50,60,70,80,90,100}/faa

# 2. Run FASconCAT for initial “all” protein concatenation
# (keeps your parallel style but now only on `faa`)
parallel -j 120 "
  cd 12-final-loci_filt/{1} && \
  perl \"$FASconCAT_DIR\"/FASconCAT-G_v1.06.1.pl -a -p -p -s -l && \
  mv FcC* ../../13-loci_concat/all/{1}/
" ::: faa

################################
# 25) THRESHOLD MATRICES (50–100%) — FAA only (36GB RAM per job)
################################
# Number of taxa (species)
export TOTAL_TAXA=$(grep -cv '^[[:space:]]*$' species.list)

# 1) Generate missing_count.txt (based on FAA presence)
find 12-final-loci_filt/faa -type f -name '*.aa.fas' -print0 \
| parallel -0 -j 120 '
  LOCUS=$(basename {} .aa.fas)
  PRES=$(grep -c "^>" 12-final-loci_filt/faa/${LOCUS}.aa.fas)
  MISSING=$(( TOTAL_TAXA - PRES ))
  printf "%s %d\n" "$LOCUS" "$MISSING"
' \
> 13-loci_concat/missing_count.txt

export TOTAL_TAXA=$(grep -cv '^[[:space:]]*$' species.list)
export FASconCAT_DIR="${FASCONCAT_DIR}"

filter_threshold() {
  local percent=$1
  local cutoff=$(( TOTAL_TAXA * (100 - percent) / 100 ))
  local matrix="13-loci_concat/matrix${percent}"
  mkdir -p "$matrix/faa"
  : > "$matrix/busco${percent}.loci.list"

  # loop over each FAA FASTA
  find 12-final-loci_filt/faa -name '*.aa.fas' -print0 \
    | while IFS= read -r -d '' file; do
        locus=$(basename "$file" .aa.fas)
        missing=$(awk -v l="$locus" '$1==l{print $2}' 13-loci_concat/missing_count.txt)
        missing=${missing:-0}
        if (( missing <= cutoff )); then
          cp -f "$file" "$matrix/faa/"
          echo "$locus" >> "$matrix/busco${percent}.loci.list"
        fi
      done

  # >>> NEW: skip if the keep-list is empty
  echo "[matrix ${percent}] no loci pass (cutoff=${cutoff}); skipping FASconCAT"
  return 0

  # Also guard against no files present (just in case)
  if compgen -G "$matrix/faa/*.aa.fas" > /dev/null; then
    ( cd "$matrix/faa" && perl "$FASconCAT_DIR"/FASconCAT-G_v1.06.1.pl -a -p -p -s -l )
  else
    echo "[matrix ${percent}] list non-empty but no files copied; skipping FASconCAT"
  fi
}
export -f filter_threshold

# run thresholds in parallel
parallel -j 120 filter_threshold ::: 50 60 70 80 90 100

##################
# 26) Summarize the Matrices — FAA only
##################
cd 13-loci_concat

# Primary (ALL)
part="all/faa/FcC_supermatrix_partition.txt"
sites=$(awk 'NF{a=$0}END{print a}' "$part" | sed -r 's/.*-(.*).*/\1/')
loci=$(wc -l < "$part")

# Threshold matrices
for percent in 50 60 70 80 90 100; do
  part="matrix${percent}/faa/FcC_supermatrix_partition.txt"
  sites=$(awk 'NF{a=$0}END{print a}' "$part" | sed -r 's/.*-(.*).*/\1/')
  loci=$(wc -l < "$part")
  echo "${percent}% matrix has ${sites} sites and ${loci} loci" | tee -a summary.extraction
done