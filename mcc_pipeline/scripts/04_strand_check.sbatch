#!/usr/bin/env bash
#SBATCH --job-name=strand_check
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/04_strand_check_%j.log
# The tutorial's flags (--rna-strandness RF / --rf / --stranded reverse) assume a TruSeq stranded
# library. Novogene's standard mRNA-seq is often unstranded, so measure it: align 2M pairs from one
# liver and one brain sample without a strand setting and count with htseq-count yes / reverse.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
W=$B/strand_check; mkdir -p "$W"; cd "$W"
N=8000000   # lines = 2M read pairs

check() {
  local s=$1
  head -n $N < <(zcat "$TRIM/${s}_1.fastq.gz") > "${s}_1.fq"   # no SIGPIPE under pipefail
  head -n $N < <(zcat "$TRIM/${s}_2.fastq.gz") > "${s}_2.fq"
  hisat2 -p "$SLURM_CPUS_PER_TASK" -x "$REF_INDEX" -1 "${s}_1.fq" -2 "${s}_2.fq" 2> "${s}.hisat2.log" |
    samtools sort -@ 4 -o "${s}.bam" -
  samtools index "${s}.bam"
  local y r
  y=$(htseq-count --format bam --order pos --mode intersection-strict --stranded yes \
        --minaqual 1 --type exon --idattr gene_id "${s}.bam" "$REF_GTF" | grep -v '^__' | awk '{s+=$2}END{print s}')
  r=$(htseq-count --format bam --order pos --mode intersection-strict --stranded reverse \
        --minaqual 1 --type exon --idattr gene_id "${s}.bam" "$REF_GTF" | grep -v '^__' | awk '{s+=$2}END{print s}')
  echo -e "$s\t$y\t$r"
}

printf "sample\tassigned_yes(FR)\tassigned_reverse(RF)\n" > strand_counts.tsv
check A1 >> strand_counts.tsv
check B1 >> strand_counts.tsv
cat strand_counts.tsv

# fraction of stranded-assigned reads that are reverse: ~1 -> RF, ~0 -> FR, ~0.5 -> unstranded
awk 'NR>1{y+=$2; r+=$3} END{f=r/(y+r); printf "reverse fraction = %.3f\n", f > "/dev/stderr";
     print (f>0.8 ? "RF" : (f<0.2 ? "FR" : "none"))}' strand_counts.tsv > "$STRAND_FILE"
echo "strandedness -> $(cat "$STRAND_FILE")"
rm -f ./*.fq
