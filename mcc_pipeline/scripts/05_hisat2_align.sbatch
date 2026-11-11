#!/usr/bin/env bash
#SBATCH --job-name=hisat2_align
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --array=0-43%22
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/05_hisat2_%A_%a.log
# rnabio.org Module 2 "Alignment": hisat2 --dta with read groups, then samtools sort.
# Piped straight into samtools sort (same result as writing the SAM first, without the 50-GB SAM).
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
S=$(sample_from_array)

# HISAT-only mode: align without --rna-strandness. HISAT2 places reads identically either way;
# the flag only sets the XS strand tag StringTie uses, and htseq-count does not read it.
HISAT_STRAND=""
# ---- STRANDEDNESS (re-enable once 04_strand_check has been run) ------------------------------
# The tutorial aligns TruSeq stranded data with --rna-strandness RF. To match it exactly, run
# 04_strand_check.sbatch, then uncomment the line below and re-align:
# strand_flags      # sets HISAT_STRAND from $STRAND_FILE (RF -> --rna-strandness RF, none -> "")
# ----------------------------------------------------------------------------------------------
mkdir -p "$ALIGN"

# PU = flowcell.lane from the first read header (@instrument:run:flowcell:lane:...)
PU=$(head -1 < <(zcat "$RAW/$S/${S}_1.fq.gz") | awk -F: '{print $3"."$4}')
TISSUE=$([[ ${S#R1_} == A* ]] && echo Liver || echo Brain)   # R1_ = Apr-2025 run

# shellcheck disable=SC2086
hisat2 -p "$SLURM_CPUS_PER_TASK" \
  --rg-id="$S" --rg "SM:$S" --rg "LB:${S}_${TISSUE}" --rg PL:ILLUMINA --rg "PU:$PU" \
  -x "$REF_INDEX" --dta $HISAT_STRAND \
  -1 "$TRIM/${S}_1.fastq.gz" -2 "$TRIM/${S}_2.fastq.gz" \
  2> "$ALIGN/${S}.hisat2.log" |
  samtools sort -@ 8 -m 2G -o "$ALIGN/${S}.bam" -
samtools index "$ALIGN/${S}.bam"
samtools flagstat -@ 4 "$ALIGN/${S}.bam" > "$ALIGN/${S}.flagstat.txt"
cat "$ALIGN/${S}.hisat2.log"
