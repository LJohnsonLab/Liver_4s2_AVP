#!/usr/bin/env python3
"""Build a custom GRCm39 reference for the EFAD / liver-APOE bulk RNA-seq.

Edits to Ensembl GRCm39 (release given on the command line):
  * mouse Apoe (ENSMUSG00000002985): removed from the GTF, gene span hard-masked (N)
  * hAPOE contig: human APOE gene region from GRCh38 (reference = E3 allele) with its
    native Ensembl exon/intron annotation. E2/E3/E4 are called afterwards from the
    BAMs at rs429358 / rs7412 rather than modelled as separate contigs.
  * Cre      : bacteriophage P1 cre CDS (GenBank X03453, protein CAA27178.1)
  * EGFP     : EGFP CDS (GenBank U55762, pEGFP-N1)
  * hAPP_5xFAD   : human APP695 CDS (NM_201414) + K670N/M671L (Swe), I716V (Fl), V717I (Lon)
  * hPSEN1_5xFAD : human PSEN1 CDS (NM_000021) + M146L, L286V
Every protein-coding insert is translated and checked before writing.
"""
import argparse, gzip, json, re, sys, urllib.request
from pathlib import Path

APOE_MOUSE = "ENSMUSG00000002985"
APOE_HUMAN = "ENSG00000130203"
PAD = 500  # bp of human flank kept on each side of APOE

CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def translate(s):
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def get(url, accept=None):
    req = urllib.request.Request(url, headers={"User-Agent": "efad-ref-build",
                                               **({"Accept": accept} if accept else {})})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode()


def ncbi_cds(acc, header_pat=None):
    """Return the (single / named) CDS nucleotide sequence of a nuccore record."""
    txt = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}"
              "&rettype=fasta_cds_na&retmode=text")
    recs = [r for r in txt.split(">") if r.strip()]
    if header_pat:
        recs = [r for r in recs if re.search(header_pat, r.splitlines()[0], re.I)]
    if len(recs) != 1:
        sys.exit(f"{acc}: expected 1 CDS, got {len(recs)}")
    return "".join(recs[0].splitlines()[1:]).upper()


def mutate(cds, pos, ref_aa, new_codon, label):
    """pos = 1-based residue in this CDS's own numbering."""
    i = (pos - 1) * 3
    old = cds[i:i + 3]
    if CODON[old] != ref_aa:
        sys.exit(f"{label}: residue {pos} is {CODON[old]} not {ref_aa}")
    return cds[:i] + new_codon + cds[i + 3:]


def check_protein(name, cds, starts_with, length):
    p = translate(cds)
    ok = p.startswith(starts_with) and p.endswith("*") and "*" not in p[:-1] and len(p) - 1 == length
    print(f"  {name}: {len(cds)} bp, {len(p)-1} aa, ORF intact={ok}")
    if not ok:
        sys.exit(f"{name} failed ORF check: {p[:30]}...{p[-10:]}")
    return p


def transgene_gtf(contig, cds, gene_name):
    L = len(cds)
    attrs = (f'gene_id "{gene_name}"; gene_version "1"; gene_name "{gene_name}"; '
             f'gene_source "custom"; gene_biotype "transgene";')
    tattrs = attrs + (f' transcript_id "{gene_name}-T1"; transcript_version "1"; '
                      f'transcript_name "{gene_name}-T1"; transcript_source "custom"; '
                      f'transcript_biotype "transgene";')
    row = lambda feat, s, e, a, frame=".": f"{contig}\tcustom\t{feat}\t{s}\t{e}\t.\t+\t{frame}\t{a}\n"
    return (row("gene", 1, L, attrs) + row("transcript", 1, L, tattrs)
            + row("exon", 1, L, tattrs + ' exon_number "1"; exon_id "' + gene_name + '-E1";')
            + row("CDS", 1, L - 3, tattrs + ' exon_number "1";', "0")
            + row("start_codon", 1, 3, tattrs + ' exon_number "1";', "0")
            + row("stop_codon", L - 2, L, tattrs + ' exon_number "1";', "0"))


def wrap(seq, n=60):
    return "\n".join(seq[i:i + n] for i in range(0, len(seq), n)) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mouse-fa", required=True)
    ap.add_argument("--mouse-gtf", required=True)
    ap.add_argument("--human-gtf", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    op = gzip.open if a.mouse_gtf.endswith(".gz") else open

    # ---- mouse Apoe span; drop it from the GTF -------------------------------------
    print("Filtering mouse GTF")
    span, kept, dropped = None, [], 0
    with op(a.mouse_gtf, "rt") as fh:
        for ln in fh:
            if not ln.startswith("#") and f'gene_id "{APOE_MOUSE}"' in ln:
                f = ln.split("\t")
                if f[2] == "gene":
                    span = (f[0], int(f[3]), int(f[4]))
                dropped += 1
                continue
            kept.append(ln)
    if not span:
        sys.exit("mouse Apoe not found")
    print(f"  mouse Apoe {span[0]}:{span[1]}-{span[2]} removed ({dropped} GTF lines)")
    overl = {re.search(r'gene_name "([^"]+)"', l).group(1) for l in kept
             if not l.startswith("#") and l.split("\t")[2] == "gene" and l.split("\t")[0] == span[0]
             and int(l.split("\t")[3]) <= span[2] and int(l.split("\t")[4]) >= span[1]
             and 'gene_name "' in l}
    print(f"  other genes overlapping masked span: {sorted(overl) or 'none'}")

    # ---- human APOE: annotation + sequence ----------------------------------------
    print("Fetching human APOE")
    hop = gzip.open if a.human_gtf.endswith(".gz") else open
    hrows = []
    with hop(a.human_gtf, "rt") as fh:
        for ln in fh:
            if not ln.startswith("#") and f'gene_id "{APOE_HUMAN}"' in ln:
                hrows.append(ln.rstrip("\n").split("\t"))
    g = next(r for r in hrows if r[2] == "gene")
    hchr, hs, he, hstrand = g[0], int(g[3]), int(g[4]), g[6]
    assert hstrand == "+", "APOE expected on + strand"
    rs, re_ = hs - PAD, he + PAD
    hseq = get(f"https://rest.ensembl.org/sequence/region/human/{hchr}:{rs}..{re_}:1"
               "?content-type=text/plain").strip().upper()
    assert len(hseq) == re_ - rs + 1
    apoe_gtf = ""
    for r in hrows:
        r = r.copy()
        r[0] = "hAPOE"; r[3] = str(int(r[3]) - rs + 1); r[4] = str(int(r[4]) - rs + 1)
        r[8] = re.sub(r'gene_name "APOE"', 'gene_name "hAPOE"', r[8])
        apoe_gtf += "\t".join(r) + "\n"
    # sanity: canonical transcript protein must carry E3 residues C112 / R158 (mature numbering)
    canon = [r for r in hrows if r[2] == "CDS" and "Ensembl_canonical" in r[8]]
    cds = "".join(hseq[int(r[3]) - rs:int(r[4]) - rs + 1] for r in sorted(canon, key=lambda r: int(r[3])))
    p = translate(cds + hseq[int(sorted(canon, key=lambda r: int(r[3]))[-1][4]) - rs + 1:][:3])
    # precursor = 18-aa signal peptide + mature; rs429358 -> residue 130, rs7412 -> 176
    print(f"  hAPOE {hchr}:{rs}-{re_} ({len(hseq)} bp); precursor aa130={p[129]} aa176={p[175]} "
          f"(E3 expects C / R)")
    assert (p[129], p[175]) == ("C", "R"), "reference APOE is not E3?"
    # positions of the two SNPs on the hAPOE contig (for genotyping later)
    def aa_to_contig(aa):  # first base of codon
        n = (aa - 1) * 3
        for r in sorted(canon, key=lambda r: int(r[3])):
            L = int(r[4]) - int(r[3]) + 1
            if n < L:
                return int(r[3]) - rs + 1 + n
            n -= L
    snps = {"rs429358 (aa130 C/R, TGC->CGC)": aa_to_contig(130),
            "rs7412 (aa176 R/C, CGC->TGC)": aa_to_contig(176)}
    for k, v in snps.items():
        print(f"  {k}: hAPOE:{v} ref={hseq[v-1]}")

    # ---- transgene CDSs -----------------------------------------------------------
    print("Fetching transgenes from NCBI")
    cre = ncbi_cds("X03453", r"protein_id=CAA27178")
    check_protein("Cre", cre, "MSNLLTVHQNLPALPVDATSDEVRKNLMD", 343)
    egfp = ncbi_cds("U55762", r"protein_id=AAB02574")
    check_protein("EGFP", egfp, "MVSKGEELFTGVVPILVELDGDVNGHKF", 239)
    app = ncbi_cds("NM_201414.2")
    check_protein("APP695 (wt)", app, "MLPGLALLLLAAWTARA", 695)
    # APP770 numbering - 75 = APP695 numbering
    app = mutate(app, 670 - 75, "K", "AAT", "APP K670N")
    app = mutate(app, 671 - 75, "M", "CTG", "APP M671L")
    app = mutate(app, 716 - 75, "I", "GTC", "APP I716V")
    app = mutate(app, 717 - 75, "V", "ATC", "APP V717I")
    pa = check_protein("hAPP_5xFAD", app, "MLPGLALLLLAAWTARA", 695)
    print(f"    Abeta region: {pa[595-1-4:642+2]}")
    ps1 = ncbi_cds("NM_000021.4")
    check_protein("PSEN1 (wt)", ps1, "MTELPAPLSYFQNAQ", 467)
    ps1 = mutate(ps1, 146, "M", "CTG", "PSEN1 M146L")
    ps1 = mutate(ps1, 286, "L", "GTG", "PSEN1 L286V")
    check_protein("hPSEN1_5xFAD", ps1, "MTELPAPLSY", 467)

    extra = {"hAPOE": hseq, "Cre": cre, "EGFP": egfp, "hAPP_5xFAD": app, "hPSEN1_5xFAD": ps1}

    # ---- write FASTA (mask Apoe) ---------------------------------------------------
    print("Writing FASTA")
    mop = gzip.open if a.mouse_fa.endswith(".gz") else open
    fa_out = out / "GRCm39_hAPOE_Cre_EGFP_5xFAD.fa"
    with mop(a.mouse_fa, "rt") as fh, open(fa_out, "w") as fo:
        name, buf = None, []

        def flush():
            if name is None:
                return
            s = "".join(buf)
            if name == span[0]:
                s = s[:span[1] - 1] + "N" * (span[2] - span[1] + 1) + s[span[2]:]
                print(f"  masked {span[2]-span[1]+1} bp on {name}")
            fo.write(f">{name}\n" + wrap(s))

        for ln in fh:
            if ln.startswith(">"):
                flush(); name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
        flush()
        for k, v in extra.items():
            fo.write(f">{k}\n" + wrap(v))
    (out / "custom_contigs.fa").write_text("".join(f">{k}\n" + wrap(v) for k, v in extra.items()))

    # ---- write GTF -----------------------------------------------------------------
    gtf_out = out / "GRCm39_hAPOE_Cre_EGFP_5xFAD.gtf"
    with open(gtf_out, "w") as fo:
        fo.writelines(kept)
        fo.write(apoe_gtf)
        for k in ("Cre", "EGFP", "hAPP_5xFAD", "hPSEN1_5xFAD"):
            fo.write(transgene_gtf(k, extra[k], k))
    (out / "apoe_snps.tsv").write_text(
        "contig\tpos\tsnp\tref\n" + "".join(f"hAPOE\t{v}\t{k.split()[0]}\t{hseq[v-1]}\n" for k, v in snps.items()))
    (out / "build_manifest.json").write_text(json.dumps({
        "mouse_apoe_masked": f"{span[0]}:{span[1]}-{span[2]}",
        "human_apoe_region": f"{hchr}:{rs}-{re_}", "apoe_snps_on_hAPOE": snps,
        "contig_lengths": {k: len(v) for k, v in extra.items()}}, indent=2))
    print(f"Done: {fa_out}\n      {gtf_out}")


if __name__ == "__main__":
    main()
