#!/usr/bin/env bash
# ============================================================================
# 01_assembly/Genome_Assembly_Single_long_CLR.sh
# ----------------------------------------------------------------------------
# Purpose:    De novo genome assembly from PacBio CLR long reads
# Inputs:     0-rawdata/*.fastq.gz
# Outputs:    Polished assemblies
# Paper ref:  Methods 'Assembly pipelines'; Suppl Table 1
# Software:   Filtlong, Flye, purge_dups, Racon, NextPolish
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

###### Genome assembly for single-end sequence, lonng read, PACBIO_SMRT Platform, CLR, PacBio RS II and Sequel II Instrument

#######Genome Assembly (WGS)


## All the raw .gz data should be in 0-rawdata folder

# I. Quality Control
## QC
mkdir -p 0-rawdata/multiqc        # for the aggregated report

# 1) Run NanoPlot in parallel, one job per core, with 4 threads each
#    Writes HTML/JSON/etc. into 0-rawdata/nanoplot/{sample}
parallel -j 15 NanoPlot \
         --fastq {} \
         --threads 4 \
         --N50 --loglength \
         -o 0-rawdata/nanoplot/{/.} \
    ::: 0-rawdata/*.fastq.gz

# 2) Run NanoStat in parallel to get tabular QC summaries
parallel -j 15 \
  'NanoStat --fastq {} > 0-rawdata/nanostat/{/.}_stats.txt' \
  ::: 0-rawdata/*.fastq.gz

# 3) (Optional) Compare multiple samples with NanoComp
FASTQS=(0-rawdata/*.fastq.gz)
TITLES=("${FASTQS[@]##*/}")
TITLES=("${TITLES[@]%.fastq.gz}")
NanoComp \
  --fastq "${FASTQS[@]}" \
  --names "${TITLES[@]}" \
  --threads 8 \
  -o 0-rawdata/nanocomp/

# 4) Aggregate all NanoPack reports with MultiQC
multiqc 0-rawdata/nanoplot/ 0-rawdata/nanostat/ 0-rawdata/nanocomp/ \
       -o 0-rawdata/multiqc/

## 3)Quality filtering with Filtlong
mkdir -p 0-1-quality_trimmed/filtlong_logs

parallel -j 7 --colsep '\s+' '
  SRAid={1}; SpeciesID={2}
  IN=0-rawdata/${SRAid}.fastq.gz
  OUT=0-1-quality_trimmed/${SRAid}.filtlong.fastq.gz
  LOG=0-1-quality_trimmed/filtlong_logs/${SpeciesID}.log

  filtlong \
    --min_length 100 \
    --min_mean_q 70 \
    --keep_percent 90 \
    $IN \
  | pigz -p 4 > $OUT 2> $LOG \
  && echo "Finished $SRAid ($SpeciesID)"
' :::: SRA_3.txt

echo "Filtlong quality filtering completed."

# 4) Assemble with Flye

parallel -j 2 --colsep '\s+' '
  SRAid={1}; SpeciesID={2}
  IN=0-1-quality_trimmed/${SRAid}.filtlong.fastq.gz
  OUT=1-assembly/${SRAid}
  LOG=1-assembly/logs/${SpeciesID}.flye.log

  flye \
    --pacbio-raw $IN \
    --threads 32 \
    --iterations 2 \
    --out-dir $OUT \
    > $LOG 2>&1 \
  && echo "Completed $SRAid ($SpeciesID)"
' :::: SRA_3.txt

echo "Flye assembly completed."

# 5) Redundancy reduction
mkdir -p 1-1-redundans

BASE=${PROJECT_ROOT}/assembly/single_long_genome_CLR

parallel -j 2 --colsep '\s+' '
  SRAid={1}

  ASM='"${BASE}"'/1-assembly/${SRAid}/assembly.fasta
  RAW='"${BASE}"'/0-rawdata/${SRAid}.fastq.gz
  OUT='"${BASE}"'/1-1-redundans/${SRAid}

  mkdir -p $OUT
  [[ -f $ASM ]] || { echo "ERROR: Missing $ASM"; exit 1; }

  (
    cd $OUT

    # 1) Map reads for coverage
    minimap2 -t8 -x map-pb "$ASM" "$RAW" | gzip -c > reads.paf.gz

    # 2) Compute histogram & cutoffs
    pbcstat reads.paf.gz
    calcuts PB.stat > cutoffs

    # 3) Self-align contigs
    split_fa "$ASM" > assembly.split.fa
    minimap2 -t8 -x asm5 -DP assembly.split.fa assembly.split.fa \
      | gzip -c > self.paf.gz

    # 4) Purge haplotigs & overlaps
    purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed

    # 5) Extract purged assembly
    get_seqs -p purged dups.bed "$ASM"
    mv purged.purged.fa assembly.purged.fasta
  )
' :::: SRA_3.txt

# 6)Polishiing

echo "=== Stage 1: Four rounds of Racon ===" | tee -a polish.log

parallel -j2 --colsep '\s+' '
  SRAid={1}

  ASSEMBLY=1-1-redundans/${SRAid}/assembly.purged.fasta
  READS=0-rawdata/${SRAid}.fastq.gz
  OUTDIR=3-polish/${SRAid}
  LOG=$OUTDIR/racon.log

  mkdir -p "$OUTDIR"
  echo "[$(date +%F\ %T)] Starting Racon for $SRAid" | tee -a "$LOG"

  # Copy initial assembly
  cp "$ASSEMBLY" "$OUTDIR/round0.fasta"

  # Four rounds of Racon
  for i in 1 2 3 4; do
    echo "[$(date +%T)] [$SRAid] Round $i: minimap2 → racon" >> "$LOG"
    minimap2 -t8 -x map-pb \
      "$OUTDIR/round$((i-1)).fasta" "$READS" \
      > "$OUTDIR/round${i}.paf" 2>> "$LOG"

    racon -t32 -m 8 -x -6 -g -8 -w 500 \
      "$READS" "$OUTDIR/round${i}.paf" \
      "$OUTDIR/round$((i-1)).fasta" \
      > "$OUTDIR/round${i}.fasta" 2>> "$LOG"
  done

  echo "[$(date +%F\ %T)] Completed Racon for $SRAid" | tee -a "$LOG" polish.log
' :::: SRA_3.txt

# Run NextPolish for final polishing
parallel -j 2 --colsep '\s+' '
  SRAid={1}

  # Define per-sample workspace
  SAMPLE_DIR=3-polish/${SRAid}
  NEXTDIR=${SAMPLE_DIR}/nextpolish_rundir
  mkdir -p ${SAMPLE_DIR}

  # 1) Prepare input list of reads
  echo "../../0-rawdata/${SRAid}.fastq.gz" > ${SAMPLE_DIR}/lgs.fofn

  # 2) Copy in your final Racon round as the draft
  cp ${SAMPLE_DIR}/round4.fasta ${SAMPLE_DIR}/draft.fasta

  # 3) Write a minimal NextPolish config
  cat > ${SAMPLE_DIR}/run.cfg <<EOF
[General]
job_type         = local
task             = best
rewrite          = yes
rerun            = 1
parallel_jobs    = 1
multithread_jobs = 32
genome           = ./draft.fasta
genome_size      = auto
workdir          = ./nextpolish_rundir
polish_options   = -p {multithread_jobs}

[lgs_option]
lgs_fofn             = ./lgs.fofn
lgs_options          = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-pb
EOF

  # 4) Run NextPolish
  cd ${SAMPLE_DIR}
  nextPolish run.cfg
' :::: SRA_3.txt

echo "NextPolish polishing completed."
