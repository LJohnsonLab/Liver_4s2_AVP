#!/usr/bin/env bash
#SBATCH --job-name=hisat2_index
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=normal
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/02_hisat2_index_%j.log
# rnabio.org Module 1 "Indexing": splice sites + exons from the GTF, then hisat2-build --ss --exon.
# (Whole-genome builds with --ss/--exon need ~160-200 GB RAM.)
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
mkdir -p "$(dirname "$REF_INDEX")"
cd "$REF"

[ -s "$ADAPTERS" ] || curl -fsSL -o "$ADAPTERS" https://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

hisat2_extract_splice_sites.py "$REF_GTF" > "$REF/splicesites.tsv"
hisat2_extract_exons.py        "$REF_GTF" > "$REF/exons.tsv"
hisat2-build -p "$SLURM_CPUS_PER_TASK" --ss "$REF/splicesites.tsv" --exon "$REF/exons.tsv" \
  "$REF_FA" "$REF_INDEX"
ls -lh "$(dirname "$REF_INDEX")"
