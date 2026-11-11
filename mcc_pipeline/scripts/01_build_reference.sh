#!/usr/bin/env bash
# Build the custom GRCm39 + hAPOE + Cre + EGFP + 5xFAD reference (FASTA + GTF).
# Runs on the MCC login node (needs internet); the STAR index itself is built by nf-core/rnaseq.
set -euo pipefail
ENS=116
B=/pscratch/lajo247_uksr/efad_bulkseq
SRC=$B/reference/src
OUT=$B/reference/GRCm39_e${ENS}_custom
mkdir -p "$SRC" "$OUT"
cd "$SRC"

get() { [ -s "$(basename "$1")" ] || curl -fsSL --retry 5 -O "$1"; }
get https://ftp.ensembl.org/pub/release-${ENS}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
get https://ftp.ensembl.org/pub/release-${ENS}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENS}.gtf.gz
get https://ftp.ensembl.org/pub/release-${ENS}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENS}.gtf.gz
get https://ftp.ensembl.org/pub/release-${ENS}/fasta/mus_musculus/dna/CHECKSUMS
sum Mus_musculus.GRCm39.dna.primary_assembly.fa.gz | awk '{print $1,$2}' > fa.sum
grep -q "$(cat fa.sum) Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" CHECKSUMS && echo "genome checksum OK"

module load python-3.9.6-gcc-8.4.1-2yf35k6 samtools-1.12-gcc-9.3.0-zo3utt7 2>/dev/null || true
python3 -I "$B/scripts/make_custom_ref.py" \
  --mouse-fa Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  --mouse-gtf Mus_musculus.GRCm39.${ENS}.gtf.gz \
  --human-gtf Homo_sapiens.GRCh38.${ENS}.gtf.gz \
  --outdir "$OUT"

samtools faidx "$OUT/GRCm39_hAPOE_Cre_EGFP_5xFAD.fa"
tail -5 "$OUT/GRCm39_hAPOE_Cre_EGFP_5xFAD.fa.fai"
cut -f1 "$OUT/GRCm39_hAPOE_Cre_EGFP_5xFAD.gtf" | grep -v '^#' | sort -u > "$OUT/gtf_seqs.txt"
comm -23 "$OUT/gtf_seqs.txt" <(cut -f1 "$OUT/GRCm39_hAPOE_Cre_EGFP_5xFAD.fa.fai" | sort) | sed 's/^/GTF contig missing from FASTA: /'
echo "reference build complete"
