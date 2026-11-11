#!/usr/bin/env bash
#SBATCH --job-name=gene_gc_len
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/pscratch/lajo247_uksr/efad_bulkseq/logs/14_gene_gc_length_%j.log
# Per-gene union-exon length and GC fraction from the custom FASTA + GTF (for cqn GC/length
# normalisation in the batch-correction benchmark).
set -euo pipefail
source /pscratch/lajo247_uksr/efad_bulkseq/scripts/config.sh
python3 -I - "$REF_FA" "$REF_GTF" "$B/results/gene_gc_length.tsv" <<'PY'
import sys, re, collections
fa, gtf, out = sys.argv[1:]
ex = collections.defaultdict(list)                     # chrom -> [(start, end, gene)]
for l in open(gtf):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    ex[f[0]].append((int(f[3]) - 1, int(f[4]), g))
res = {}
def process(chrom, seq):
    by_gene = collections.defaultdict(list)
    for s, e, g in ex.get(chrom, []): by_gene[g].append((s, e))
    for g, iv in by_gene.items():
        iv.sort(); merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
            else: merged.append([s, e])
        L = gc = 0
        for s, e in merged:
            sub = seq[s:e].upper(); L += len(sub); gc += sub.count("G") + sub.count("C")
        res[g] = (L, gc / L if L else float("nan"))
name, buf = None, []
for l in open(fa):
    if l.startswith(">"):
        if name: process(name, "".join(buf))
        name, buf = l[1:].split()[0], []
    else: buf.append(l.strip())
if name: process(name, "".join(buf))
with open(out, "w") as fo:
    fo.write("GeneID\tlength\tgc\n")
    fo.writelines(f"{g}\t{L}\t{gc:.4f}\n" for g, (L, gc) in res.items())
print(len(res), "genes")
PY
