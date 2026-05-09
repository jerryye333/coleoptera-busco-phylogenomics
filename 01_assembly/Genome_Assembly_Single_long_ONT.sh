#!/usr/bin/env bash
# ============================================================================
# 01_assembly/Genome_Assembly_Single_long_ONT.sh
# ----------------------------------------------------------------------------
# Purpose:    De novo genome assembly from Oxford Nanopore long reads
# Inputs:     0-rawdata/*.fastq.gz
# Outputs:    Polished assemblies
# Paper ref:  Methods 'Assembly pipelines'; Suppl Table 1
# Software:   Porechop, Filtlong, Flye, purge_dups, Racon, Medaka
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

###### Genome assembly for single-end sequins, long read, OXFORD_NANOPORE Platform, ONT, MinION and PromethION Instrument

#######Genome Assembly (WGS)

## All the raw .gz data should be in 0-rawdata folder
# Quality Control
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

# I. Adapter trimming

# Create output directory
OUTDIR=0-1-adapter_trimmed
mkdir -p "${OUTDIR}"

# Parallel over SRA_4.txt: col1=SRAid, col2=SpeciesID
parallel -j 2 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  INFILE=0-rawdata/${SRAid}.fastq.gz
  OUTFILE='"${OUTDIR}"'/${SRAid}_trimmed.fastq.gz

  # Run porechop if input exists

  porechop \
    --input "${INFILE}" \
    --output "${OUTFILE}" \
    --threads 32 \
    --verbosity 2 \
    2>&1
' :::: SRA_4.txt

# II. Quality filtering with Filtlong
mkdir -p 0-2-quality_trimmed/filtlong_logs

parallel -j 5 --colsep '\s+' '
  SRAid={1}; SpeciesID={2}
  IN=0-1-adapter_trimmed/${SRAid}_trimmed.fastq.gz
  OUT=0-2-quality_trimmed/${SRAid}.filtlong.fastq.gz
  LOG=0-2-quality_trimmed/filtlong_logs/${SpeciesID}.log

  filtlong \
    --min_length 100 \
    --min_mean_q 70 \
    --keep_percent 90 \
    $IN \
  | pigz -p 4 > $OUT 2> $LOG \
  && echo "Finished $SRAid ($SpeciesID)"
' :::: SRA_4.txt

echo "Filtlong quality filtering completed."

# III. Genome Assembly with Flye

parallel -j 2 --colsep '\s+' '
  SRAid={1}; SpeciesID={2}
  IN=0-2-quality_trimmed/${SRAid}.filtlong.fastq.gz
  OUT=1-assembly/${SRAid}
  LOG=1-assembly/logs/${SpeciesID}.flye.log

  flye \
    --nano-raw $IN \
    --threads 32 \
    --iterations 2 \
    --scaffold \
    --out-dir $OUT \
    > $LOG 2>&1 \
  && echo "Completed $SRAid ($SpeciesID)"
' :::: SRA_4.txt

echo "Flye assembly completed."

# IV. Redundancy reduction
mkdir -p 1-1-redundans

BASE=${PROJECT_ROOT}/assembly/single_long_genome_ONT

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
    minimap2 -t8 -x map-ont "$ASM" "$RAW" | gzip -c > reads.paf.gz

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
' :::: SRA_4.txt

# V) Polishiing

echo "=== Stage 1: Four rounds of Racon ===" | tee -a polish.log

parallel -j2 --colsep '\s+' '
  SRAid={1}

  ASSEMBLY=1-1-redundans/${SRAid}/assembly.purged.fasta
  READS=0-rawdata/${SRAid}.fastq.gz
  OUTDIR=2-polish/${SRAid}
  LOG=$OUTDIR/racon.log

  mkdir -p "$OUTDIR"
  echo "[$(date +%F\ %T)] Starting Racon for $SRAid" | tee -a "$LOG"

  # Copy initial assembly
  cp "$ASSEMBLY" "$OUTDIR/round0.fasta"

  # Four rounds of Racon
  for i in 1 2 3 4; do
    echo "[$(date +%T)] [$SRAid] Round $i: minimap2 → racon" >> "$LOG"
    minimap2 -t8 -x map-ont \
      "$OUTDIR/round$((i-1)).fasta" "$READS" \
      > "$OUTDIR/round${i}.paf" 2>> "$LOG"

    racon -t32 -m 8 -x -6 -g -8 -w 500 \
      "$READS" "$OUTDIR/round${i}.paf" \
      "$OUTDIR/round$((i-1)).fasta" \
      > "$OUTDIR/round${i}.fasta" 2>> "$LOG"
  done

  echo "[$(date +%F\ %T)] Completed Racon for $SRAid" | tee -a "$LOG" polish.log
' :::: SRA_4.txt

# VI) Medaka polishing
echo "=== Stage 2: Medaka polishing ===" | tee -a polish.log

parallel -j4 --colsep '\s+' '
  SRAid={1}

  READS=0-rawdata/${SRAid}.fastq.gz
  ROUND4=2-polish/${SRAid}/round4.fasta
  OUTDIR=2-polish/${SRAid}
  LOG=$OUTDIR/medaka.log

  echo "[$(date +%F\ %T)] Starting Medaka for $SRAid" | tee -a "$LOG"

  medaka_consensus \
    -i "$READS" \
    -d "$ROUND4" \
    -o "$OUTDIR/medaka_out" \
    -t 32 \
    -m r941_min_high_g507 \
    2>> "$LOG"

  mv "$OUTDIR/medaka_out/consensus.fasta" "$OUTDIR/polished.fasta"
  echo "[$(date +%F\ %T)] Completed Medaka for $SRAid" | tee -a "$LOG" polish.log
' :::: SRA_4.txt
