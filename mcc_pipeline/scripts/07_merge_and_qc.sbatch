#!/usr/bin/env bash
#SBATCH --job-name=merge_qc
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/07_merge_qc_%j.log
# Merge per-sample outputs into matrices (as in rnabio.org Module 3) and run MultiQC.
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
OUT=$B/results; mkdir -p "$OUT"
mapfile -t S < "$SAMPLES"

# ---- htseq-count matrix: gene_read_counts_table_all_final.tsv (drop __ summary rows) ----
# Same table the tutorial builds with `join`, merged in Python so every sample is checked to have
# the identical gene rows (htseq appends __ rows after the genes, which breaks join's sort order).
python3 -I - "$EXPR/htseq_counts" "$OUT" "${S[@]}" <<'EOF'
import sys
src, out, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
cols = {}
for s in samples:
    with open(f"{src}/{s}_gene.tsv") as fh:
        cols[s] = dict(l.rstrip("\n").split("\t") for l in fh)
keys = list(cols[samples[0]])
assert all(list(cols[s]) == keys for s in samples), "htseq gene rows differ between samples"
genes = [k for k in keys if not k.startswith("__")]
summ  = [k for k in keys if k.startswith("__")]
def write(path, head, rows):
    with open(path, "w") as fo:
        fo.write(head + "\t" + "\t".join(samples) + "\n")
        fo.writelines(f"{k}\t" + "\t".join(cols[s][k] for s in samples) + "\n" for k in rows)
write(f"{out}/gene_read_counts_table_all_final.tsv", "GeneID", genes)
write(f"{out}/htseq_summary.tsv", "category", summ)   # no_feature / ambiguous / too_low_aQual / ...
print(f"count matrix: {len(genes)} genes x {len(samples)} samples")
EOF

# ---- gene ID -> name / biotype from the custom GTF (tutorial's ENSG_ID2Name.txt) ----
awk -F'\t' '$3=="gene"{id=""; nm=""; bt="";
  if (match($9,/gene_id "[^"]+"/))      id=substr($9,RSTART+9,RLENGTH-10);
  if (match($9,/gene_name "[^"]+"/))    nm=substr($9,RSTART+11,RLENGTH-12);
  if (match($9,/gene_biotype "[^"]+"/)) bt=substr($9,RSTART+14,RLENGTH-15);
  if (nm=="") nm=id; print id"\t"nm"\t"bt"\t"$1}' "$REF_GTF" > "$OUT/gene_id2name.tsv"

# ---- StringTie gene TPM / FPKM matrices from gene_abundances.tsv ----
python3 -I - "$EXPR/stringtie/ref_only" "$OUT" "${S[@]}" <<'EOF'
import sys, csv
src, out, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
for metric in ("TPM", "FPKM"):
    mat, names = {}, {}
    for s in samples:
        with open(f"{src}/{s}/gene_abundances.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gid = r["Gene ID"]; names[gid] = r["Gene Name"]
                mat.setdefault(gid, {})[s] = mat.get(gid, {}).get(s, 0.0) + float(r[metric])
    with open(f"{out}/gene_{metric.lower()}_all_samples.tsv", "w") as fo:
        fo.write("GeneID\tGeneName\t" + "\t".join(samples) + "\n")
        for g in sorted(mat):
            fo.write(f"{g}\t{names[g]}\t" + "\t".join(f"{mat[g].get(s, 0):.4f}" for s in samples) + "\n")
print("stringtie matrices written")
EOF

# ---- per-sample QC table: fastp + hisat2 + htseq-count summary -> sample_qc_metrics.tsv ----
python3 -I - "$QC/fastp" "$ALIGN" "$OUT/htseq_summary.tsv" "$OUT" "${S[@]}" <<'EOF'
import sys, json, re, csv
fp, al, ht, out, samples = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
htrows = {r["category"]: r for r in csv.DictReader(open(ht), delimiter="\t")}
counts = {}
for line in open(f"{out}/gene_read_counts_table_all_final.tsv"):
    if line.startswith("GeneID"): hdr = line.rstrip("\n").split("\t")[1:]; continue
    for s, v in zip(hdr, line.rstrip("\n").split("\t")[1:]): counts[s] = counts.get(s, 0) + int(v)
def pct(pat, txt):
    m = re.search(pat, txt); return float(m.group(1)) if m else float("nan")
rows = []
for s in samples:
    j = json.load(open(f"{fp}/{s}.fastp.json"))
    before, after = j["summary"]["before_filtering"], j["summary"]["after_filtering"]
    h = open(f"{al}/{s}.hisat2.log").read()
    hs = {k.lstrip("_"): int(htrows[k][s]) for k in htrows}
    assigned = counts[s]
    total_ht = assigned + sum(hs.values())
    rows.append({
        "sample": s,
        "raw_read_pairs": before["total_reads"] // 2,
        "trimmed_read_pairs": after["total_reads"] // 2,
        "pct_reads_passing_fastp": round(100 * after["total_reads"] / before["total_reads"], 2),
        "q30_rate_after": round(100 * after["q30_rate"], 2),
        "gc_after": round(100 * after["gc_content"], 2),
        "duplication_rate": round(100 * j.get("duplication", {}).get("rate", float("nan")), 2),
        "adapter_trimmed_reads": j.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
        "insert_size_peak": j.get("insert_size", {}).get("peak", float("nan")),
        "hisat2_overall_alignment_pct": pct(r"([\d.]+)% overall alignment rate", h),
        "hisat2_concordant_unique_pct": pct(r"\(([\d.]+)%\) aligned concordantly exactly 1 time", h),
        "hisat2_concordant_multi_pct": pct(r"\(([\d.]+)%\) aligned concordantly >1 times", h),
        "htseq_assigned": assigned,
        **{f"htseq_{k}": v for k, v in hs.items()},
        "pct_assigned_to_genes": round(100 * assigned / total_ht, 2) if total_ht else float("nan"),
    })
with open(f"{out}/sample_qc_metrics.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"sample QC table: {len(rows)} samples")
EOF

# ---- MultiQC over FastQC (raw + trimmed), fastp, hisat2, samtools, htseq ----
multiqc -f -o "$QC/multiqc" -n efad_bulkseq_multiqc \
  "$QC/fastqc_raw" "$QC/fastqc_trimmed" "$QC/fastp" "$ALIGN" "$EXPR/htseq_counts"
ls -lh "$OUT"
