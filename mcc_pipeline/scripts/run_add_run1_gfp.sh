#!/usr/bin/env bash
# Interim: process the 20 April-2025 GFP samples (R1_, already on MCC) now and merge with the
# 44 current-run samples (64 total). The 2025-04 Cre samples are excluded (bad data).
set -euo pipefail
B=/pscratch/lajo247_uksr/efad_bulkseq
GFP_IDS="A1 A2 A3 A4 A5 A11 A12 A13 A14 A15 B1 B2 B3 B4 B5 B11 B12 B13 B14 B15"
{ cat "$B/samples.txt"; printf 'R1_%s\n' $GFP_IDS; } > "$B/samples_gfp.txt"
export SAMPLES=$B/samples_gfp.txt
source "$B/scripts/config.sh"
cd "$B/scripts"
[ "$(wc -l < "$SAMPLES")" -eq 64 ] || { echo "samples_gfp.txt must list 64 samples"; exit 1; }
sub() { sbatch --parsable --export=ALL,SAMPLES="$SAMPLES" "$@"; }
FP=$(sub --array=44-63%20 03_fastp.sbatch)
AL=$(sub --array=44-63%20 --dependency=afterok:$FP 05_hisat2_align.sbatch)
EX=$(sub --array=44-63%20 --dependency=afterok:$AL 06_expression.sbatch)
MQ=$(sub --dependency=afterok:$EX 07_merge_and_qc.sbatch)
AP=$(sub --dependency=afterok:$MQ 08_apoe_genotype_transgenes.sbatch)
echo "$EX" > "$LOGS/run1_gfp_expression_jobid.txt"
echo "run1 GFP: fastp $FP | align $AL | expression $EX | merge(64) $MQ | apoe $AP" | tee "$LOGS/run_add_run1_gfp_jobids.txt"
