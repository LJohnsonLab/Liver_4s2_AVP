#!/usr/bin/env bash
#SBATCH --job-name=apoe_geno
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/08_apoe_genotype_%j.log
# Sample-identity checks (not part of the tutorial):
#  (1) APOE rs429358 / rs7412 base counts on the hAPOE contig -> E2/E3/E4 call per sample
#      (the same idea as rnabio.org "Alignment Visualization - Allele Expression")
#  (2) Cre / EGFP / hAPP / hPSEN1 / hAPOE rows from the htseq-count matrix.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
O=$B/results/apoe_transgene_check; mkdir -p "$O"
mapfile -t SAMP < "$SAMPLES"

printf "sample\tsnp\tdepth\tT\tC\n" > "$O/apoe_snp_counts.tsv"
tail -n +2 "$REF/apoe_snps.tsv" | while IFS=$'\t' read -r ctg pos snp ref; do
  printf '%s\n' "${SAMP[@]}" | while read -r s; do
    # samtools mpileup base string -> strip read-start (^ + MAPQ char), read-end ($) and indels
    samtools mpileup -r "$ctg:$pos-$pos" -f "$REF_FA" -Q 20 -q 20 -d 1000000 "$ALIGN/${s}.bam" 2>/dev/null |
      python3 -I -c '
import sys, re
s, snp = sys.argv[1], sys.argv[2]
line = sys.stdin.readline().split("\t")
if len(line) < 5: print(f"{s}\t{snp}\t0\t0\t0"); sys.exit()
ref, bases = line[2].upper(), line[4]
bases = re.sub(r"\^.", "", bases).replace("$", "")
while (m := re.search(r"[+-](\d+)", bases)):           # drop inserted/deleted sequence
    n = int(m.group(1)); bases = bases[:m.start()] + bases[m.end() + n:]
b = bases.upper().replace(".", ref).replace(",", ref)
n = sum(b.count(x) for x in "ACGT"); t = b.count("T"); c = b.count("C")
print(f"{s}\t{snp}\t{n}\t{t}\t{c}")' "$s" "$snp"
  done
done >> "$O/apoe_snp_counts.tsv"

# rs429358: T = Cys (E2/E3), C = Arg (E4);  rs7412: C = Arg (E3/E4), T = Cys (E2)
python3 -I - "$O" <<'EOF'
import sys, csv, collections
o = sys.argv[1]
d = collections.defaultdict(dict)
for r in csv.DictReader(open(f"{o}/apoe_snp_counts.tsv"), delimiter="\t"):
    t, c = int(r["T"]), int(r["C"])
    d[r["sample"]][r["snp"]] = (t, c)
def call(e4, n4, e2, n2):
    if n4 < 10 or n2 < 10: return "low_coverage"
    a4 = "E4/E4" if e4 > .8 else ("E3/E4" if e4 > .2 else None)
    a2 = "E2/E2" if e2 > .8 else ("E2/E3" if e2 > .2 else None)
    if a4 and a2: return "E2/E4"
    return a4 or a2 or "E3/E3"
with open(f"{o}/apoe_genotype_calls.tsv", "w") as fo:
    fo.write("sample\trs429358_C_frac\trs429358_n\trs7412_T_frac\trs7412_n\tcall\n")
    for s in sorted(d, key=lambda x: (x.startswith("R1_"), x.removeprefix("R1_")[0], int(x.removeprefix("R1_")[1:]))):
        t4, c4 = d[s].get("rs429358", (0, 0)); t2, c2 = d[s].get("rs7412", (0, 0))
        n4, n2 = t4 + c4, t2 + c2
        e4 = c4 / n4 if n4 else 0; e2 = t2 / n2 if n2 else 0
        fo.write(f"{s}\t{e4:.2f}\t{n4}\t{e2:.2f}\t{n2}\t{call(e4, n4, e2, n2)}\n")
print(open(f"{o}/apoe_genotype_calls.tsv").read())
EOF

M=$B/results/gene_read_counts_table_all_final.tsv
{ head -1 "$M"; grep -P '^(Cre|EGFP|hAPP_5xFAD|hPSEN1_5xFAD|ENSG00000130203)\t' "$M"; } > "$O/construct_counts.tsv"
cat "$O/construct_counts.tsv"
