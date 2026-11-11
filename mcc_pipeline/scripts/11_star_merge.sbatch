#!/usr/bin/env bash
#SBATCH --job-name=star_merge
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/11_star_merge_%j.log
# Merge STAR ReadsPerGene (column 2 = unstranded, matching htseq-count --stranded no) and
# Log.final.out mapping stats into results/.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
mapfile -t S < "$SAMPLES"
python3 -I - "$STAR_DIR" "$B/results" "${S[@]}" <<'PY'
import sys
src, out, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
counts, summ, keys = {}, {}, None
for s in samples:
    rows = [l.rstrip("\n").split("\t") for l in open(f"{src}/{s}/ReadsPerGene.out.tab")]
    counts[s] = {r[0]: r[1] for r in rows if not r[0].startswith("N_")}
    summ[s]   = {r[0]: r[1] for r in rows if r[0].startswith("N_")}
    k = [r[0] for r in rows if not r[0].startswith("N_")]
    assert keys is None or k == keys, f"gene rows differ in {s}"; keys = k
with open(f"{out}/star_gene_counts.tsv", "w") as fo:
    fo.write("GeneID\t" + "\t".join(samples) + "\n")
    fo.writelines(g + "\t" + "\t".join(counts[s][g] for s in samples) + "\n" for g in keys)
fields = ["Number of input reads", "Uniquely mapped reads %", "% of reads mapped to multiple loci",
          "% of reads unmapped: too short", "% of reads unmapped: other", "Mismatch rate per base, %"]
with open(f"{out}/star_mapping_summary.tsv", "w") as fo:
    fo.write("sample\t" + "\t".join(f.replace(" ", "_").replace("%", "pct").replace(",", "").replace(":", "")
                                    for f in fields) + "\tN_noFeature\tN_ambiguous\tN_multimapping\n")
    for s in samples:
        log = dict(l.strip().split(" |\t") for l in open(f"{src}/{s}/Log.final.out") if " |\t" in l)
        fo.write(s + "\t" + "\t".join(log.get(f, "NA").rstrip("%") for f in fields) + "\t" +
                 "\t".join(summ[s][k] for k in ("N_noFeature", "N_ambiguous", "N_multimapping")) + "\n")
print(f"STAR counts: {len(keys)} genes x {len(samples)} samples")
PY
