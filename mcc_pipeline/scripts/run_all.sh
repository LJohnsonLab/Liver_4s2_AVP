#!/usr/bin/env bash
# Submit the whole pipeline with SLURM dependencies. Run on the MCC login node AFTER
# 00_verify_md5.sh passes and 01_build_reference.sh has written the custom FASTA/GTF.
#   index ─┐
#   fastp ─┴─> align (array) ─> expression (array) ─> merge + MultiQC ─> APOE / construct check
# HISAT-only mode: the strand check is not in the chain; everything runs unstranded (see config.sh).
# Skips the index / fastp / alignment when their outputs already exist.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
cd "$B/scripts"
test -s "$REF_FA" && test -s "$REF_GTF" || { echo "reference not built"; exit 1; }
[ "$(wc -l < "$SAMPLES")" -eq 44 ] || { echo "samples.txt must list 44 samples"; exit 1; }

sub() { sbatch --parsable "$@"; }
# index: reuse if built, wait on it if already queued/running, else submit
if [ -s "${REF_INDEX}.1.ht2" ]; then IDX=""
elif J=$(squeue -u "$USER" -h -n hisat2_index -o %i | head -1) && [ -n "$J" ]; then IDX=$J
else IDX=$(sub 02_hisat2_index.sbatch); fi
# fastp: skip if every sample already has both trimmed reads
if [ "$(ls "$TRIM"/*_[12].fastq.gz 2>/dev/null | wc -l)" -eq 88 ]; then FP=""; else FP=$(sub 03_fastp.sbatch); fi
DEP=$(echo "$IDX $FP" | xargs | tr ' ' ':')

# ---- STRANDEDNESS CHECK (disabled; re-enable to match the tutorial's stranded settings) -------
# Uncomment to measure strandedness before alignment, and swap the AL line below for the
# commented one so alignment waits on it. 06/07/08 then use the measured setting automatically.
# SC=$(sub ${DEP:+--dependency=afterok:$DEP} 04_strand_check.sbatch)
# DEP=$SC
# -----------------------------------------------------------------------------------------------

# alignment: skip if every sample already has an indexed BAM
if [ "$(ls "$ALIGN"/*.bam.bai 2>/dev/null | wc -l)" -eq 44 ]; then AL=""
else AL=$(sub ${DEP:+--dependency=afterok:$DEP} 05_hisat2_align.sbatch); fi
EX=$(sub ${AL:+--dependency=afterok:$AL} 06_expression.sbatch)
MQ=$(sub --dependency=afterok:$EX 07_merge_and_qc.sbatch)
AP=$(sub --dependency=afterok:$MQ 08_apoe_genotype_transgenes.sbatch)
printf "index %s | fastp %s | align %s | expression %s | merge %s | apoe %s | strandedness %s\n" \
  "$IDX" "$FP" "$AL" "$EX" "$MQ" "$AP" "$( [ -s "$STRAND_FILE" ] && cat "$STRAND_FILE" || echo none)" \
  | tee "$LOGS/run_all_jobids.txt"
