#!/usr/bin/env bash
#SBATCH --job-name=star_align
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/10_star_%A_%a.log
# STAR on the same fastp-trimmed reads used for HISAT2. Only gene counts are kept
# (--quantMode GeneCounts, --outSAMtype None): this is an internal comparison, no BAMs needed.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
S=$(sample_from_array)
O=$STAR_DIR/$S; mkdir -p "$O"
STAR --runThreadN "$SLURM_CPUS_PER_TASK" --genomeDir "$STAR_INDEX" \
  --readFilesIn "$TRIM/${S}_1.fastq.gz" "$TRIM/${S}_2.fastq.gz" --readFilesCommand zcat \
  --outSAMtype None --quantMode GeneCounts \
  --outTmpDir "$O/_tmp" --outFileNamePrefix "$O/"
cat "$O/Log.final.out"
