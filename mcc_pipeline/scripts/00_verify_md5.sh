#!/usr/bin/env bash
#SBATCH --job-name=md5_check
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/00_md5_%j.log
# Verify transferred FASTQs against Novogene's per-sample MD5.txt (16 samples in parallel).
set -uo pipefail
F=/pscratch/lajo247_uksr/efad_bulkseq/fastq
cd "$F"
ls -d */ | tr -d / | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -I{} \
  sh -c 'cd "{}" && md5sum -c --quiet MD5.txt >/dev/null 2>&1 && echo "OK {}" || echo "FAIL {}"' \
  | sort -V > "$F/../logs/md5_results.txt"
n=$(wc -l < "$F/../logs/md5_results.txt"); bad=$(grep -c '^FAIL' "$F/../logs/md5_results.txt")
grep '^FAIL' "$F/../logs/md5_results.txt"
echo "checked $n sample folders; failures: $bad"
[ "$bad" -eq 0 ] && [ "$n" -eq 44 ]
