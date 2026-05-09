#!/usr/bin/env bash
# ============================================================================
# 01_assembly/Genome_Assembly_Pair-end_short.sh
# ----------------------------------------------------------------------------
# Purpose:    De novo genome assembly from paired-end short reads (Illumina)
# Inputs:     0-rawdata/*.fastq.gz; SRA_1.txt (SRA_id<TAB>species_id)
# Outputs:    3-gapclosing/<species>_scaffolds.gapcloser.fa
# Paper ref:  Methods 'Assembly pipelines'; Suppl Table 1
# Software:   fastp, BBTools, MEGAHIT, Redundans, BESST, GapCloser
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

#######Genome Assembly (WGS) paired-end short read


# For redundans, you need to create a new conda environment with python 3.10

# For BESST, you need to create a new conda environment

cd pair-end_short_assembly
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
' :::: SRA_1.txt

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

' :::: SRA_1.txt

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
' :::: SRA_1.txt

echo "clumpify.sh finished."

mkdir -p 0-3-duplicate_removed/fastqc
mkdir -p 0-3-duplicate_removed/multiqc
parallel -j 20 fastqc -t 4 {} -o 0-3-duplicate_removed/fastqc ::: 0-3-duplicate_removed/*.fq.gz
multiqc 0-3-duplicate_removed/fastqc/ -o 0-3-duplicate_removed/multiqc/

# II. Genome Assembly
## multi-kmer genome assembly using megahit
mkdir -p 1-assembly
touch 1-assembly/assembly.statistics

# Run MEGAHIT per species in parallel
parallel -j 4 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  R1=0-3-duplicate_removed/${SRAid}_1.clumped.fq.gz
  R2=0-3-duplicate_removed/${SRAid}_2.clumped.fq.gz
  OUTDIR=1-assembly/${SpeciesID}_megahit

  megahit -1 $R1 -2 $R2 \
          -o $OUTDIR \
          --k-list 21,41,61,81,101,121 \
          -t 32 \
          -m 0.25 \
          > 1-assembly/logs/${SpeciesID}.log 2>&1

  if [[ -f "$OUTDIR/final.contigs.fa" ]]; then
    statswrapper.sh in=$OUTDIR/final.contigs.fa format=6 >> 1-assembly/assembly.statistics
  else
    echo "No assembly output for $SpeciesID ($SRAid)" >> log.txt
  fi
' :::: SRA_1.txt

## Reduction of heterozygous contigs using redundans

parallel -j 8 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN=1-assembly/${SpeciesID}_megahit/final.contigs.fa
  OUTDIR=1-assembly/${SpeciesID}_megahit/reduced
  REDUCED_FINAL=1-assembly/${SpeciesID}_megahit/final.contigs.reduced.fa
  LOG=1-assembly/redundans_logs/${SpeciesID}.redundans.log

  redundans.py -v -f $IN -o $OUTDIR -t 16 \
               --log $LOG \
               --noscaffolding --nogapclosing \
               --identity 0.85 \
               --minimap2reduce

  if [[ -f "$OUTDIR/scaffolds.reduced.fa" ]]; then
    mv $OUTDIR/scaffolds.reduced.fa $REDUCED_FINAL
  else
    echo "Redundans failed for $SpeciesID ($SRAid): no output file found."
  fi

  rm -rf $OUTDIR *fa.fai log* reads.list
' :::: SRA_1.txt

echo "Redundans reduction step finished."

# III. Scaffolding
## Scaffolding assembled contigs
echo "Scaffolding......"
mkdir -p 2-scaffolding

parallel -j 2 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  CONTIGS=1-assembly/${SpeciesID}_megahit/final.contigs.reduced.fa
  R1=0-3-duplicate_removed/${SRAid}_1.clumped.fq.gz
  R2=0-3-duplicate_removed/${SRAid}_2.clumped.fq.gz
  OUTBAM=2-scaffolding/${SpeciesID}.mapped.bam
  THREADS=16
  LOGFILE=2-scaffolding/bam_logs/${SpeciesID}.log

  minimap2 -ax sr "$CONTIGS" "$R1" "$R2" -t $THREADS \
    | samtools view -bS -@ $THREADS \
    | samtools sort -@ $THREADS -m 2G -o "$OUTBAM"

  samtools index "$OUTBAM" -@ $THREADS
' :::: SRA_1.txt

## Scaffolding using BESST

# Run BESST scaffolding with default parameters
mkdir -p 2-scaffolding/BESST_output

parallel -j 4 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  # Paths
  CONTIGS=1-assembly/${SpeciesID}_megahit/final.contigs.reduced.fa
  BAM=2-scaffolding/${SpeciesID}.mapped.bam
  OUTDIR=2-scaffolding/BESST_output/${SpeciesID}
  SCAFF="$OUTDIR/BESST_output/pass1/Scaffolds_pass1.fa"
  LOGFILE=2-scaffolding/BESST_logs/${SpeciesID}.log

  runBESST -c "$CONTIGS" -f "$BAM" -o "$OUTDIR" -orientation fr --iter 10000 > "$LOGFILE" 2>&1

  if [[ ! -s "$SCAFF" ]]; then
    echo "BESST scaffolding FAILED for $SRAid ($SpeciesID)"
    exit 1
  fi
' :::: SRA_1.txt

echo "BESST scaffolding step finished."

# IV. Gap Filling
## Gap filling

mkdir -p 3-gapclosing
mkdir -p 3-gapclosing/bb_logs

parallel -j 30 gunzip -k {} ::: 0-3-duplicate_removed/*.fq.gz

parallel -j 4 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  # Generate a per-species config file for GapCloser
  LIBCONFIG=3-gapclosing/gapcloser_config_{2}.cfg
  echo "[LIB]" > $LIBCONFIG
  echo "q1=0-3-duplicate_removed/${SRAid}_1.clumped.fq" >> $LIBCONFIG
  echo "q2=0-3-duplicate_removed/${SRAid}_2.clumped.fq" >> $LIBCONFIG

  # Define input scaffold file from BESST output (using the SpeciesID from {2})
  IN=2-scaffolding/BESST_output/{2}/BESST_output/pass1/Scaffolds_pass1.fa
  # Define output file for gap closing
  OUT=3-gapclosing/{2}_scaffolds.gapcloser.fa
  THREADS=10
  READ_LENGTH=155

  GapCloser -a "$IN" -b $LIBCONFIG -o "$OUT" -t 16 -l $READ_LENGTH -t $THREADS
  if [[ ! -s "$OUT" ]]; then
    exit 1
  fi
' :::: SRA_1.txt
