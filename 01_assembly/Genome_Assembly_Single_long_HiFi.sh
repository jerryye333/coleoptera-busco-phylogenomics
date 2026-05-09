#!/usr/bin/env bash
# ============================================================================
# 01_assembly/Genome_Assembly_Single_long_HiFi.sh
# ----------------------------------------------------------------------------
# Purpose:    De novo genome assembly from PacBio HiFi long reads
# Inputs:     0-rawdata/*.fastq.gz; SRA_2.txt
# Outputs:    2-1-redundans/genome/genome_<species>.fa
# Paper ref:  Methods 'Assembly pipelines'; Suppl Table 1
# Software:   SeqKit, Hifiasm, purge_dups
# ----------------------------------------------------------------------------
# Set PROJECT_ROOT before running any section.
# ============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

###### Genome assembly for single-end sequence, lonng read, PACBIO_SMRT Platform, HiFi, Revio and Sequel IIe Instrument

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

# ## Quality filtering using seqkit
mkdir -p 1-seqkit_filtered

parallel -j 5 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}

  IN="0-rawdata/${SRAid}.fastq.gz"
  OUT="1-seqkit_filtered/${SRAid}.filt.fastq.gz"
  LOG="1-seqkit_filtered/seqkit_logs/${SpeciesID}.log"


  seqkit seq "$IN" \
    -j 16 \
    --min-qual 20 \
    -o "$OUT" \
    2> "$LOG"

  if [[ ! -s "$OUT" ]]; then
    echo "ERROR: Empty output for $SRAid ($SpeciesID)"
  fi
' :::: SRA_2.txt

echo "SeqKit quality filtering finished."

# II. Genome Assembly
## HiFi assembly with hifiasm
mkdir -p 2-assembly/hifiasm

# Run up to 4 assemblies in parallel (4 jobs × 32 threads = 128 cores)
parallel -j 2 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}
  IN=1-seqkit_filtered/${SRAid}.filt.fastq.gz
  OUTDIR=2-assembly/hifiasm/${SpeciesID}
  LOG=2-assembly/logs/${SpeciesID}_hifiasm.log
  GFA=$OUTDIR/assembly.bp.p_ctg.gfa

# Can use --dual-scaf for scaffolding
  mkdir -p "$OUTDIR"
  hifiasm -t32 \
  -i \
  -f0 \
  -D 10.0 \
  -N 500 \
  --min-hist-cnt 10 \
  -o "$OUTDIR/assembly" \
  "$IN" > "$LOG" 2>&1

  # Convert GFA to FASTA
  awk '"'"'/^S/ {print ">" $2; print $3}'"'"' \
    "$GFA" > "$OUTDIR/assembly.p_ctg.fa"
' :::: SRA_2.txt

# Assemble for large genomes with hifiasm
parallel -j 1 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}
  IN=1-seqkit_filtered/${SRAid}.filt.fastq.gz
  OUTDIR=2-assembly/hifiasm/${SpeciesID}
  LOG=2-assembly/logs/${SpeciesID}_hifiasm.log
  GFA=$OUTDIR/assembly.bp.p_ctg.gfa

  mkdir -p "$OUTDIR"
  hifiasm -t32 \
  -i \
  -f38 \
  -D 5.0 \
  -N 100 \
  -r 1 \
  --hg-size 3.26g \
  --min-hist-cnt 10 \
  -o "$OUTDIR/assembly" \
  "$IN" > "$LOG" 2>&1

  # Convert GFA to FASTA
  awk '"'"'/^S/ {print ">" $2; print $3}'"'"' \
    "$GFA" > "$OUTDIR/assembly.p_ctg.fa"
' :::: SRA_2_remain.txt

# III. Redundancy reduction (optional; apply when duplicated BUSCOs are high)
mkdir -p 2-1-redundans

BASE=${PROJECT_ROOT}/assembly/single_long_genome_HiFi

parallel -j 2 --colsep '\s+' '
  SRAid={1}
  SpeciesID={2}
  ASM='"${BASE}"'/2-assembly/hifiasm/${SpeciesID}/genome_${SpeciesID}.fa
  RAW='"${BASE}"'/0-rawdata/${SRAid}.fastq.gz
  OUT='"${BASE}"'/2-1-redundans/${SRAid}

  mkdir -p $OUT
  [[ -f $ASM ]] || { echo "ERROR: Missing $ASM"; exit 1; }

  (
    cd $OUT

    # 1) Map reads for coverage
    minimap2 -t8 -x map-hifi "$ASM" "$RAW" | gzip -c > reads.paf.gz

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
' :::: SRA_2.txt
