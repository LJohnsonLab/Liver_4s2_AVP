#!/usr/bin/env bash
# STAR internal comparison for every analysed sample (samples_gfp.txt: 44 current-run + 20
# April-2025 GFP). Waits for the R1 GFP fastp jobs if they are still queued.
set -euo pipefail
B=/pscratch/lajo247_uksr/efad_bulkseq
export SAMPLES=$B/samples_gfp.txt
source "$B/scripts/config.sh"
cd "$B/scripts"
N=$(wc -l < "$SAMPLES")
sub() { sbatch --parsable --export=ALL,SAMPLES="$SAMPLES" "$@"; }
if [ -s "$STAR_INDEX/SA" ]; then IX=""; else IX=$(sub 09_star_index.sbatch); fi
AL=$(sub --array=0-$((N - 1))%16 ${IX:+--dependency=afterok:$IX} 10_star_align.sbatch)
MG=$(sub --dependency=afterok:$AL 11_star_merge.sbatch)
echo "STAR: index ${IX:-reused} | align $AL | merge $MG" | tee "$LOGS/run_star_jobids.txt"
