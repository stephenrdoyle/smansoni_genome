"""Z-W gametologue divergence in Trichobilharzia regenti (WBPS19, PRJEB44434).

Genes on candidate W contigs (unplaced contigs absent from the treg_1 male,
smansoni_genome.z_strata_remap.md section 4) are paired with their closest
gene on CHR_Z, codon-aligned and given dN/dS as in zw_ds.py. Each pair is
also placed on S. mansoni, and where S. mansoni has its own Z-W pair for the
same gene, a four-point test on protein distances asks whether the two W
copies group together (recombination stopped before the lineages split) or
each W groups with its own Z (stopped independently).

usage: treg_zw.py TREG_DATA SMAN_DATA SMAN_ZW_TSV OUT_DIR [W_CONTIG ...]
"""
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, _count_diff_NG86, _count_site_NG86, _get_codon_list, cal_dn_ds

TABLE = CodonTable.generic_by_id[1]
treg_dir, sman_dir, sman_zw, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
W_CONTIGS = set(sys.argv[5:]) or {"UNPLACED_00115", "UNPLACED_00117", "UNPLACED_00373", "UNPLACED_00423", "UNPLACED_00478"}
out.mkdir(exist_ok=True)


def read_genes(path):
    genes = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            genes[attrs["ID"].removeprefix("gene:")] = (f[0], int(f[3]), attrs.get("biotype", ""))
    return genes


def longest_with_cds(prot_path, cds_path, gene_of):
    best = {}
    with gzip.open(prot_path, "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            g = gene_of(r)
            if g not in best or len(r.seq) > len(best[g].seq):
                best[g] = r
    with gzip.open(cds_path, "rt") as fh:
        cds_all = {r.id: r for r in SeqIO.parse(fh, "fasta")}
    prot = {g: str(r.seq).rstrip("*") for g, r in best.items() if r.id in cds_all}
    cds = {g: str(cds_all[best[g].id].seq).upper() for g in prot}
    return prot, cds


def gene_attr(r):
    """WBPS fasta headers carry gene=<id>; fall back to the transcript id stem."""
    for tok in r.description.split():
        if tok.startswith("gene="):
            return tok.removeprefix("gene=")
    return r.id.rsplit(".", 1)[0]


# --- helpers as in zw_ds.py -----------------------------------------------------
def kmers(s, k):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def build_index(seqs, k):
    idx = defaultdict(set)
    for g, s in seqs.items():
        for km in kmers(s, k):
            idx[km].add(g)
    return idx


def top_candidates(q, index, n=15):
    hits = Counter()
    for km in q:
        for g in index.get(km, ()):
            hits[g] += 1
    return [g for g, _ in hits.most_common(n)]


blosum = substitution_matrices.load("BLOSUM62")
local_aa = PairwiseAligner(mode="local", substitution_matrix=blosum, open_gap_score=-11, extend_gap_score=-1)
global_aa = PairwiseAligner(mode="global", substitution_matrix=blosum, open_gap_score=-11, extend_gap_score=-1,
                            end_gap_score=0)


def aligned_pairs(aln):
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            yield a0 + k, b0 + k


def codon_strings(w_nt, z_nt, codon_map):
    wc, zc = [], []
    for zi in sorted(codon_map):
        w, z = codon_map[zi], z_nt[3 * zi:3 * zi + 3]
        if len(z) < 3 or len(w) < 3 or "N" in w or "N" in z:
            continue
        if str(Seq(w).translate()) == "*" or str(Seq(z).translate()) == "*":
            continue
        wc.append(w)
        zc.append(z)
    return "".join(wc), "".join(zc)


def estimate(w_str, z_str):
    res = {"codons": len(w_str) // 3}
    if res["codons"] < 30:
        return res
    a, b = CodonSeq(w_str), CodonSeq(z_str)
    for method in ("NG86", "YN00", "ML"):
        try:
            dn, ds = cal_dn_ds(a, b, method=method)
        except Exception:
            dn, ds = float("nan"), float("nan")
        res[f"dN_{method}"], res[f"dS_{method}"] = dn, ds
    l1, l2 = _get_codon_list(a), _get_codon_list(b)
    S = (_count_site_NG86(l1, TABLE)[0] + _count_site_NG86(l2, TABLE)[0]) / 2
    N = (_count_site_NG86(l1, TABLE)[1] + _count_site_NG86(l2, TABLE)[1]) / 2
    sd = nd = 0
    for x, y in zip(l1, l2):
        d = _count_diff_NG86(x, y, TABLE)
        sd, nd = sd + d[0], nd + d[1]
    res["pS"], res["pN"] = sd / S, nd / N
    res["nt_identity"] = sum(x == y for x, y in zip(w_str, z_str)) / len(w_str)
    return res


def aa_dist(p, q):
    """p-distance over positions aligned in a global protein alignment."""
    pairs = list(aligned_pairs(global_aa.align(p, q)[0]))
    return 1 - sum(p[i] == q[j] for i, j in pairs) / max(len(pairs), 1), len(pairs)


# --- data -----------------------------------------------------------------------------
tpre = "trichobilharzia_regenti"
spre = "schistosoma_mansoni.PRJEA36577.WBPS19"
tgenes = read_genes(treg_dir / f"{tpre}.annotations.gff3.gz")
tprot, tcds = longest_with_cds(treg_dir / f"{tpre}.protein.fa.gz", treg_dir / f"{tpre}.CDS_transcripts.fa.gz", gene_attr)
sgenes = read_genes(sman_dir / f"{spre}.annotations.gff3.gz")
sprot, _ = longest_with_cds(sman_dir / f"{spre}.protein.fa.gz", sman_dir / f"{spre}.CDS_transcripts.fa.gz",
                            lambda r: r.id.split(".")[0])
# S. mansoni Z-W pairs (full-length rows) from zw_ds.py
szw = pd.read_csv(sman_zw, sep="\t")
sman_w_of_z = dict(zip(szw.loc[szw.region == "full", "Z_gene"], szw.loc[szw.region == "full", "W_gene"]))

w_genes = sorted(g for g, v in tgenes.items() if v[0] in W_CONTIGS and g in tprot)
print(f"{len(w_genes)} annotated genes on {len(W_CONTIGS)} candidate W contigs", file=sys.stderr)
others = {g: s for g, s in tprot.items() if tgenes.get(g, ("",))[0] not in W_CONTIGS}
t_idx = build_index(others, 4)
s_idx = build_index(sprot, 4)

rows = []
for w in w_genes:
    wp, w_nt = tprot[w], tcds[w]
    scored = sorted(((local_aa.score(wp, others[g]), g) for g in top_candidates(kmers(wp, 4), t_idx)), reverse=True)
    if not scored:
        continue
    on_z = [(sc, g) for sc, g in scored if tgenes[g][0] == "CHR_Z"]
    best_sc, best_g = scored[0]
    z_sc, z = on_z[0] if on_z else (0, "")
    row = {"W_gene": w, "W_contig": tgenes[w][0], "W_aa": len(wp), "best_hit": best_g,
           "best_hit_chrom": tgenes[best_g][0], "best_score": best_sc, "Z_gene": z, "Z_score": z_sc,
           "Z_is_best": z == best_g}
    if z:
        zp, z_nt = tprot[z], tcds[z]
        pairs = list(aligned_pairs(global_aa.align(wp, zp)[0]))
        row["aa_identity"] = sum(wp[i] == zp[j] for i, j in pairs) / max(len(pairs), 1)
        row["W_aa_coverage"], row["Z_aa_coverage"] = len(pairs) / len(wp), len(pairs) / len(zp)
        row.update(estimate(*codon_strings(w_nt, z_nt, {j: w_nt[3 * i:3 * i + 3] for i, j in pairs})))
        # S. mansoni homologue of the Z copy, and its own W copy if it has one
        ss = sorted(((local_aa.score(zp, sprot[g]), g) for g in top_candidates(kmers(zp, 4), s_idx)), reverse=True)
        if ss:
            sm = ss[0][1]
            row.update({"Sm_gene": sm, "Sm_chrom": sgenes.get(sm, ("",))[0], "Sm_start": sgenes.get(sm, ("", 0))[1]})
            sw = sman_w_of_z.get(sm)
            if sw and sw in sprot:
                d = {k: aa_dist(*pq)[0] for k, pq in {
                    "TrW_TrZ": (wp, zp), "SmW_SmZ": (sprot[sw], sprot[sm]), "TrW_SmW": (wp, sprot[sw]),
                    "TrZ_SmZ": (zp, sprot[sm]), "TrW_SmZ": (wp, sprot[sm]), "TrZ_SmW": (zp, sprot[sw])}.items()}
                split = {"independent (TrW,TrZ | SmW,SmZ)": d["TrW_TrZ"] + d["SmW_SmZ"],
                         "shared W (TrW,SmW | TrZ,SmZ)": d["TrW_SmW"] + d["TrZ_SmZ"],
                         "crossed (TrW,SmZ | TrZ,SmW)": d["TrW_SmZ"] + d["TrZ_SmW"]}
                row.update({"Sm_W_gene": sw, **{f"d_{k}": v for k, v in d.items()},
                            "four_point": min(split, key=split.get)})
    rows.append(row)

df = pd.DataFrame(rows)
if "dS_ML" in df:
    df["dNdS_ML"] = df["dN_ML"] / df["dS_ML"]
df.to_csv(out / "treg_zw_gametologue_ds.tsv", sep="\t", index=False, float_format="%.4g")
cols = [c for c in ["W_gene", "W_contig", "Z_gene", "Z_is_best", "best_hit_chrom", "aa_identity", "codons",
                    "pS", "dS_ML", "dN_ML", "dNdS_ML", "Sm_gene", "Sm_chrom", "Sm_start", "Sm_W_gene", "four_point"]
        if c in df]
print(df[cols].to_string(index=False))
