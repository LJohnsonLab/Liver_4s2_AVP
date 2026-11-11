# Shared settings for the EFAD / liver-APOE bulk RNA-seq pipeline (Griffith lab rnabio.org workflow)
B=/pscratch/lajo247_uksr/efad_bulkseq
ENV=$B/env/rnaseq
export PATH=$ENV/bin:$PATH

REF=$B/reference/GRCm39_e116_custom
REF_FA=$REF/GRCm39_hAPOE_Cre_EGFP_5xFAD.fa
REF_GTF=$REF/GRCm39_hAPOE_Cre_EGFP_5xFAD.gtf
REF_INDEX=$REF/hisat2/GRCm39_custom          # hisat2 index prefix
ADAPTERS=$B/reference/illumina_multiplex.fa  # genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

RAW=$B/fastq                 # <sample>/<sample>_{1,2}.fq.gz (Novogene layout)
TRIM=$B/trimmed
QC=$B/qc
ALIGN=$B/alignments/hisat2
EXPR=$B/expression
LOGS=$B/logs
SAMPLES=${SAMPLES:-$B/samples.txt}   # one ID per line; samples_all.txt adds the R1_ (Apr-2025 run) GFP samples
STRAND_FILE=$B/strandedness.txt   # written by 04_strand_check: RF | FR | none

sample_from_array() { sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES"; }

# Library strandedness -> per-tool flags (tutorial uses RF / --rf / reverse for TruSeq stranded)
# HISAT-only mode: strandedness is not measured yet, so everything runs unstranded ("none").
# Unstranded counting is valid for either library type (it only loses the ability to separate
# overlapping antisense genes). Once 04_strand_check has written $STRAND_FILE, rerun
# 06_expression -> 07 -> 08 and the measured setting is picked up automatically.
strand_flags() {
  local s="none"
  [ -s "$STRAND_FILE" ] && s=$(cat "$STRAND_FILE")
  case $s in
    RF)   HISAT_STRAND="--rna-strandness RF"; STRINGTIE_STRAND="--rf"; HTSEQ_STRAND="reverse" ;;
    FR)   HISAT_STRAND="--rna-strandness FR"; STRINGTIE_STRAND="--fr"; HTSEQ_STRAND="yes" ;;
    none) HISAT_STRAND="";                    STRINGTIE_STRAND="";     HTSEQ_STRAND="no" ;;
    *) echo "bad strandedness '$s' in $STRAND_FILE" >&2; return 1 ;;
  esac
}

# STAR (internal comparison with the HISAT2 / htseq-count counts; not part of the tutorial)
STAR_INDEX=$REF/star          # built by 09_star_index.sbatch from the same custom FASTA + GTF
STAR_DIR=$B/alignments/star
