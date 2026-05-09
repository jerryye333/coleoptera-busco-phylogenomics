#!/usr/bin/env bash
# ============================================================================
# 01_assembly/Transcriptome_Assembly.sh
# ----------------------------------------------------------------------------
# Purpose:    De novo transcriptome assembly from paired-end Illumina RNA-seq
# Inputs:     0-rawdata/*.fastq.gz
# Outputs:    Clustered transcript assemblies
# Paper ref:  Methods 'Assembly pipelines'; Suppl Table 1
# Software:   fastp, BBTools, rnaSPAdes, CD-HIT
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

########### Transriptome Assembly

# For redundans, you need to create a new conda environment with python 3.10


# I. Quality Control
## fastqc
mkdir -p 0-rawdata/fastqc
mkdir -p 0-rawdata/multiqc

parallel -j 10 fastqc {} -o 0-rawdata/fastqc ::: 0-rawdata/*.fastq.gz
echo "Fastqc for rawdata finished."

# Generate MultiQC report for raw data
multiqc 0-rawdata/fastqc/ -o 0-rawdata/multiqc/

## adapter trimming

mkdir -p 0-1-adapter_trimmed
parallel -j 8 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN1=0-rawdata/${SRAid}_1.fastq.gz
  IN2=0-rawdata/${SRAid}_2.fastq.gz
  OUT1=0-1-adapter_trimmed/${SRAid}_1_adapttrimmed.fq.gz
  OUT2=0-1-adapter_trimmed/${SRAid}_2_adapttrimmed.fq.gz

  fastp \
    -i "$IN1" -I "$IN2" \
    -o "$OUT1" -O "$OUT2" \
    -Q -L --detect_adapter_for_pe \
    --thread 16 \
    -j 0-1-adapter_trimmed/fastp_reports/${SpeciesID}.json \
    -h 0-1-adapter_trimmed/fastp_reports/${SpeciesID}.html \
    > 0-1-adapter_trimmed/fastp_reports/${SpeciesID}.log 2>&1
' :::: SRA_6.txt

## quality trimming

parallel -j 8 --colsep '\s+' '
  SRAid={1}; SpeciesID={2}

  IN1=0-1-adapter_trimmed/${SRAid}_1_adapttrimmed.fq.gz
  IN2=0-1-adapter_trimmed/${SRAid}_2_adapttrimmed.fq.gz
  TMP1=0-2-quality_trimmed/${SRAid}_1.tmp.fq.gz
  TMP2=0-2-quality_trimmed/${SRAid}_2.tmp.fq.gz
  OUT1=0-2-quality_trimmed/${SRAid}_1.quatrim.fq.gz
  OUT2=0-2-quality_trimmed/${SRAid}_2.quatrim.fq.gz
  LOG=0-2-quality_trimmed/bbduk_logs/${SpeciesID}.log

  # 1) Trim + correct
  bbduk.sh \
    in1=$IN1 in2=$IN2 \
    out1=$TMP1 out2=$TMP2 \
    pigz=16 ordered qtrim=rl trimq=20 minlen=15 \
    ecco=t maxns=5 \
    trimpolya=10 trimpolyg=10 trimpolyc=10 \
    threads=16 \
    > $LOG 2>&1

  # 2) Repair pairs in-place
  repair.sh \
    in1=$TMP1 in2=$TMP2 \
    out1=$OUT1 out2=$OUT2 \
    overwrite=t \
    outs=0-2-quality_trimmed/repair_stats/${SpeciesID}.stats \
    repair

  # 3) Clean up temps if desired
  rm $TMP1 $TMP2

' :::: SRA_6.txt

## remove duplicates
mkdir -p 0-3-duplicate_removed
parallel -j 5 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN1=0-2-quality_trimmed/${SRAid}_1.quatrim.fq.gz
  IN2=0-2-quality_trimmed/${SRAid}_2.quatrim.fq.gz
  OUT1=0-3-duplicate_removed/${SRAid}_1.clumped.fq.gz
  OUT2=0-3-duplicate_removed/${SRAid}_2.clumped.fq.gz
  LOG=0-3-duplicate_removed/clumpify_logs/${SpeciesID}.log

  echo "Clumpifying $SRAid ($SpeciesID)..."

  clumpify.sh in1=$IN1 in2=$IN2 out1=$OUT1 out2=$OUT2 \
              deletetemp=t \
              -Xmx20g threads=16 pigz=16 dedupe \
              > $LOG 2>&1
' :::: SRA_6.txt

echo "clumpify.sh finished."

mkdir -p 0-3-duplicate_removed/fastqc
mkdir -p 0-3-duplicate_removed/multiqc
parallel -j 20 fastqc -t 4 {} -o 0-3-duplicate_removed/fastqc ::: 0-3-duplicate_removed/*.fq.gz
multiqc 0-3-duplicate_removed/fastqc/ -o 0-3-duplicate_removed/multiqc/

# II. Transcriptome Assembly (using rnaSPAdes)
mkdir 1-assembly
parallel -j 5 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  # Define the input read files for the sample
  R1=0-3-duplicate_removed/${SRAid}_1.clumped.fq.gz
  R2=0-3-duplicate_removed/${SRAid}_2.clumped.fq.gz

  # Define the output directory for rnaspades assembly
  OUTDIR=1-assembly/{2}

    rnaspades.py -1 "$R1" -2 "$R2" -o "$OUTDIR" -t 20 -m 64
' :::: SRA_6.txt

echo "rnaSPAdes assembly step finished."

## Reduction of heterozygous contigs using redundans

parallel -j 8 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN=1-assembly/${SpeciesID}/transcripts.fasta
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
' :::: SRA_6.txt

echo "Redundans reduction step finished."

# Or use cd-hit-est to reduce redundancy (another round if the duplicate is still high, recommended)
mkdir -p output

parallel -j 8 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN=${PROJECT_ROOT}/assembly/assembled_transcriptome/transcriptome_${SpeciesID}.fasta
  OUTDIR=output/${SpeciesID}/cdhit90
  NR_FINAL=output/${SpeciesID}/transcriptome_${SpeciesID}.fasta
  LOG=output/cdhit_logs/${SpeciesID}.cdhit90.log

  mkdir -p "$(dirname "$LOG")" "$OUTDIR"

  cd-hit-est \
    -i "$IN" \
    -o "$OUTDIR/cluster90" \
    -c 0.90 \
    -n 8 \
    -aS 0.90 \
    -G 0 \
    -g 1 \
    -d 0 \
    -T 16 \
    -M 16000 \
    > "$LOG" 2>&1

  if [[ -f "$OUTDIR/cluster90" ]]; then
    mv "$OUTDIR/cluster90" "$NR_FINAL"
    # move the cluster report too (useful for QC)
    if [[ -f "$OUTDIR/cluster90.clstr" ]]; then
      mv "$OUTDIR/cluster90.clstr" "${NR_FINAL%.fa}.clstr"
    fi
    echo "Non-redundant transcripts saved to $NR_FINAL"
  else
    echo "CD-HIT-EST failed for $SpeciesID ($SRAid): no output file found."
  fi
' :::: species.txt
