#!/usr/bin/env bash
#SBATCH --job-name=icre_merge
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/13_icre_merge_%j.log
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
mapfile -t S < "$SAMPLES"
{ printf "sample\tcontig\treads\n"; for s in "${S[@]}"; do cat "$B/expression/construct_check/${s}.tsv"; done; } \
  > "$B/results/apoe_transgene_check/construct_unaligned_reads.tsv"
awk -F'\t' 'NR>1 && $2=="iCre_AY056050"' "$B/results/apoe_transgene_check/construct_unaligned_reads.tsv"
