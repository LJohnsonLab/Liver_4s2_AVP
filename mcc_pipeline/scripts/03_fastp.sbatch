#!/usr/bin/env bash
#SBATCH --job-name=fastp
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-43%22
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/03_fastp_%A_%a.log
# rnabio.org Module 1 "Pre-alignment QC" (FastQC on raw reads) +
# Module 2 "Adapter Trim" (fastp: illumina_multiplex.fa, 13-bp hard left trim, min length 25).
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
S=$(sample_from_array)
mkdir -p "$TRIM" "$QC/fastqc_raw/$S" "$QC/fastqc_trimmed/$S" "$QC/fastp"

fastqc -t 2 --outdir "$QC/fastqc_raw/$S" "$RAW/$S/${S}_1.fq.gz" "$RAW/$S/${S}_2.fq.gz"

fastp -i "$RAW/$S/${S}_1.fq.gz" -I "$RAW/$S/${S}_2.fq.gz" \
      -o "$TRIM/${S}_1.fastq.gz" -O "$TRIM/${S}_2.fastq.gz" \
      -l 25 --adapter_fasta "$ADAPTERS" --trim_front1 13 --trim_front2 13 \
      --thread "$SLURM_CPUS_PER_TASK" \
      --json "$QC/fastp/${S}.fastp.json" --html "$QC/fastp/${S}.fastp.html" \
      2> "$QC/fastp/${S}.fastp.log"

fastqc -t 2 --outdir "$QC/fastqc_trimmed/$S" "$TRIM/${S}_1.fastq.gz" "$TRIM/${S}_2.fastq.gz"
echo "done $S"
