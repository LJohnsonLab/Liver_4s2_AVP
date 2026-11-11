#!/usr/bin/env bash
#SBATCH --job-name=star_index
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/09_star_index_%j.log
# STAR genome index from the SAME custom FASTA + GTF as the HISAT2 index. Reads are 137 bp after
# fastp (150 - 13), so sjdbOverhang = 136.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
mkdir -p "$STAR_INDEX"
STAR --runMode genomeGenerate --runThreadN "$SLURM_CPUS_PER_TASK" \
  --genomeDir "$STAR_INDEX" --genomeFastaFiles "$REF_FA" \
  --sjdbGTFfile "$REF_GTF" --sjdbOverhang 136 \
  --outTmpDir "$STAR_INDEX/_tmp" --outFileNamePrefix "$STAR_INDEX/"
