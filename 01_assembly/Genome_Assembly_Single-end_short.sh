#!/usr/bin/env bash
# ============================================================================
# 01_assembly/Genome_Assembly_Single-end_short.sh
# ----------------------------------------------------------------------------
# Purpose:    De novo genome assembly from single-end short reads (Illumina NextSeq 500)
# Inputs:     0-rawdata/*.fastq.gz
# Outputs:    Scaffolded assemblies (gap-filled)
# Paper ref:  Methods 'Assembly pipelines'; Suppl Table 1
# Software:   fastp, BBTools, MEGAHIT, Redundans, ABySS-Sealer
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

###### Genome assembly for single-end sequences, short read, ILLUMINA Platform, NextSeq 500 Instrument

#######Genome Assembly (WGS)

# For redundans, you need to create a new conda environment with python 3.10

## All the raw .gz data should be in 0-rawdata folder

# I. Quality Control
## fastqc
mkdir -p 0-rawdata/fastqc
mkdir -p 0-rawdata/multiqc

# Run FastQC on each compressed file in 0-rawdata, writing output to 0-rawdata/fastqc/
parallel -j 20 fastqc -t 4 {} -o 0-rawdata/fastqc ::: 0-rawdata/*.fastq.gz

# Generate MultiQC report for raw data
multiqc 0-rawdata/fastqc/ -o 0-rawdata/multiqc/

## adapter trimming

mkdir -p 0-1-adapter_trimmed

parallel -j 8 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN=0-rawdata/${SRAid}.fastq.gz
  OUT=0-1-adapter_trimmed/${SRAid}_adapttrimmed.fq.gz

  fastp \
    -i "$IN" \
    -o "$OUT" \
    -Q -L \
    --thread 16 \
    --trim_poly_g \
    -j 0-1-adapter_trimmed/fastp_reports/${SpeciesID}.json \
    -h 0-1-adapter_trimmed/fastp_reports/${SpeciesID}.html \
    > 0-1-adapter_trimmed/fastp_reports/${SpeciesID}.log 2>&1
' :::: SRA_5.txt

## quality trimming

parallel -j 8 --colsep '\s+' '
  SRAid={1}; SpeciesID={2}

  IN=0-1-adapter_trimmed/${SRAid}_adapttrimmed.fq.gz
  OUT=0-2-quality_trimmed/${SRAid}.quatrim.fq.gz
  LOG=0-2-quality_trimmed/bbduk_logs/${SpeciesID}.log

  bbduk.sh \
    in="$IN" \
    out="$OUT" \
    pigz=16 ordered \
    qtrim=rl trimq=20 minlen=15 \
    ecco=t maxns=5 \
    trimpolya=10 trimpolyg=10 trimpolyc=10 \
    threads=16 \
    int=f \
    overwrite=t \
    > "$LOG" 2>&1
' :::: SRA_5.txt

## remove duplicates
mkdir -p 0-3-duplicate_removed
parallel -j 5 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN=0-2-quality_trimmed/${SRAid}.quatrim.fq.gz
  OUT=0-3-duplicate_removed/${SRAid}.clumped.fq.gz
  LOG=0-3-duplicate_removed/clumpify_logs/${SpeciesID}.log

  echo "Clumpifying $SRAid ($SpeciesID)..."

 clumpify.sh \
    in="$IN" \
    out="$OUT" \
    deletetemp=t \
    -Xmx20g \
    threads=16 \
    pigz=16 \
    dedupe \
    > "$LOG" 2>&1
' :::: SRA_5.txt

echo "clumpify.sh finished."

mkdir -p 0-3-duplicate_removed/fastqc
mkdir -p 0-3-duplicate_removed/multiqc
parallel -j 20 fastqc -t 4 {} -o 0-3-duplicate_removed/fastqc ::: 0-3-duplicate_removed/*.fq.gz
multiqc 0-3-duplicate_removed/fastqc/ -o 0-3-duplicate_removed/multiqc/

# II. Genome Assembly
## genome assembly using SPAdes
mkdir -p 1-assembly

parallel -j 5 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  # Single‐end trimmed reads
  R=0-3-duplicate_removed/${SRAid}.clumped.fq.gz
  OUTDIR=1-assembly/${SpeciesID}

  spades.py \
  --isolate \
    -s "$R" \
    -o "$OUTDIR" \
    -t 20 \
    -m 64 \
    > 1-assembly/logs/${SpeciesID}.spades.log 2>&1
' :::: SRA_5.txt

## Reduction of heterozygous contigs using redundans

parallel -j 8 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN=1-assembly/${SpeciesID}/scaffolds.fasta
  OUTDIR=1-assembly/${SpeciesID}/reduced
  REDUCED_FINAL=1-assembly/${SpeciesID}/final.reduced.fa
  LOG=1-assembly/redundans_logs/${SpeciesID}.redundans.log

  redundans.py -v -f $IN -o $OUTDIR -t 16 \
               --log $LOG \
               --noscaffolding --nogapclosing \
               --identity 0.85

  if [[ -f "$OUTDIR/scaffolds.reduced.fa" ]]; then
    mv $OUTDIR/scaffolds.reduced.fa $REDUCED_FINAL
  else
    echo "Redundans failed for $SpeciesID ($SRAid): no output file found."
  fi

  rm -rf $OUTDIR *fa.fai log* reads.list
' :::: SRA_5.txt

echo "Redundans reduction step finished."

# III. Gap Filling

mkdir -p 2-gapfilling

parallel -j 5 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  ASSEMBLY=1-assembly/${SpeciesID}/final.reduced.fa
  READS=0-3-duplicate_removed/${SRAid}.clumped.fq.gz
  OUT=2-gapfilling_single-end/${SpeciesID}_contig.sealed
  BLOOM_SIZE="20G"  # Adjust based on genome size
  THREADS=20

    abyss-sealer \
      -k21 -k41 -k61 -k81 -k101 \
      -b "$BLOOM_SIZE" \
      -t "$THREADS" \
      -o "$OUT" \
      -S "$ASSEMBLY" \
      "$READS"
    if [[ -s "${OUT}_scaffold.fa" ]]; then
        echo "Success: ${SpeciesID}"
    else
        exit 1
    fi
' :::: SRA_5.txt

echo "Sealer gap filling step finished."

