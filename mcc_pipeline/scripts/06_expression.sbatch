#!/usr/bin/env bash
#SBATCH --job-name=expression
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=0-43%22
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/06_expression_%A_%a.log
# rnabio.org Module 3 "Expression": StringTie reference-only (-e -B) for TPM/FPKM and
# htseq-count (intersection-strict, minaqual 1, exon, gene_id) for raw counts -> DESeq2.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
strand_flags
S=$(sample_from_array)
mkdir -p "$EXPR/stringtie/ref_only/$S" "$EXPR/htseq_counts"

# shellcheck disable=SC2086
stringtie $STRINGTIE_STRAND -p "$SLURM_CPUS_PER_TASK" -G "$REF_GTF" -e -B \
  -o "$EXPR/stringtie/ref_only/$S/transcripts.gtf" \
  -A "$EXPR/stringtie/ref_only/$S/gene_abundances.tsv" \
  "$ALIGN/${S}.bam"

htseq-count --format bam --order pos --mode intersection-strict --stranded "$HTSEQ_STRAND" \
  --minaqual 1 --type exon --idattr gene_id \
  "$ALIGN/${S}.bam" "$REF_GTF" > "$EXPR/htseq_counts/${S}_gene.tsv"
tail -5 "$EXPR/htseq_counts/${S}_gene.tsv"
