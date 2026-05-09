#!/usr/bin/env bash
# ============================================================================
# 02_busco/run_busco.sh
# ----------------------------------------------------------------------------
# Purpose:    Run BUSCO v6 (coleoptera_odb12), then gather single+fragment orthologs for extraction
# Inputs:     Genome / transcriptome assemblies; species lists
# Outputs:    busco_result/run_<species>/run_coleoptera_odb12/busco_sequences/*.faa
# Paper ref:  Methods 'BUSCO assessment'; Fig 1; Suppl Fig 1
# Software:   BUSCO v6.0.0, GNU parallel
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

# ============================================================================
# 1) Run BUSCO v6 against coleoptera_odb12 (3,729 loci)
# ============================================================================
export _JAVA_OPTIONS="-Xms1g -Xmx8g"

# Genome mode
parallel -j 4 '
  species="{}"
  fa="${PROJECT_ROOT}/assemblies/genome_${species}.fa"
  conda run -n busco busco -i "$fa" \
    -l ${PROJECT_ROOT}/busco_downloads/lineages/coleoptera_odb12 \
    -m geno -c 30 \
    --out_path ${PROJECT_ROOT}/BUSCO/genome \
    -o "run_${species}"
' :::: species_genome.txt

# Transcriptome mode (same database)
parallel -j 4 '
  species="{}"
  fa="${PROJECT_ROOT}/assemblies/transcriptome_${species}.fa"
  conda run -n busco busco -i "$fa" \
    -l ${PROJECT_ROOT}/busco_downloads/lineages/coleoptera_odb12 \
    -m transcriptome -c 30 \
    --out_path ${PROJECT_ROOT}/BUSCO/transcriptome \
    -o "run_${species}"
' :::: species_transcriptome.txt

# ============================================================================
# 2) Gather single + fragmented BUSCO sequences into busco_result/<species>/
#    (Paper uses single-copy + fragmented orthologs; duplicates excluded.)
# ============================================================================
parallel -j 30 '
  species={}
  base="${PROJECT_ROOT}/BUSCO/genome/run_${species}/run_coleoptera_odb12/busco_sequences"
  dest="${PROJECT_ROOT}/busco_result/run_${species}/run_coleoptera_odb12/busco_sequences"
  mkdir -p "$dest"
  for sub in single_copy_busco_sequences fragmented_busco_sequences; do
    src="$base/$sub"
    [[ -d "$src" ]] && cp -a --update=none --reflink=auto "$src"/. "$dest"/
  done
' :::: species.txt
