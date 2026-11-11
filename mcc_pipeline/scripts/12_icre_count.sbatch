#!/usr/bin/env bash
#SBATCH --job-name=icre_count
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/12_icre_%A_%a.log
# The AAV is EGFP-T2A-iCre. iCre (codon-improved Cre, GenBank AY056050) differs from the
# wild-type P1 Cre contig in the genome reference, so iCre reads end up unaligned. Re-align each
# sample's unaligned reads to a small index (iCre + EGFP + wild-type Cre) and count them.
# Gene counts are unaffected (these reads never reached a mouse gene).
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
CK=$B/reference/construct_check
S=$(sample_from_array)
OUT=$B/expression/construct_check; mkdir -p "$OUT"
if [ ! -s "$CK/construct.1.ht2" ]; then      # first task builds the index; others wait for it
  if mkdir "$CK/.lock" 2>/dev/null; then
    samtools faidx "$REF/custom_contigs.fa"
    { cat "$CK/iCre.fa"; samtools faidx "$REF/custom_contigs.fa" EGFP Cre; } > "$CK/construct.fa"
    hisat2-build -q "$CK/construct.fa" "$CK/construct"
  else
    until [ -s "$CK/construct.8.ht2" ]; do sleep 10; done; sleep 5
  fi
fi
# unaligned primary reads (both mates) -> single-end, unspliced alignment to the construct index
samtools fastq -@ 4 -f 4 -F 0x900 "$ALIGN/${S}.bam" 2>/dev/null |
  hisat2 -p "$SLURM_CPUS_PER_TASK" -x "$CK/construct" -U - --no-spliced-alignment --no-unal \
    2> "$OUT/${S}.hisat2.log" |
  samtools view -q 1 - | cut -f3 | sort | uniq -c |
  awk -v s="$S" '{print s"\t"$2"\t"$1}' > "$OUT/${S}.tsv"
# total primary reads for scaling
awk -v s="$S" 'NR==2{print s"\t__primary_reads\t"$1}' "$ALIGN/${S}.flagstat.txt" >> "$OUT/${S}.tsv"
cat "$OUT/${S}.tsv"
