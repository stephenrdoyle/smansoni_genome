"""Z-W gametologue divergence (dN, dS) for S. mansoni v10 (WBPS19).

Pairs every gene in SM_V10_WSR with its best match elsewhere in the
genome, codon-aligns the pair and estimates dN/dS. Protein-coding W genes
are paired by protein similarity; W pseudogenes by nucleotide similarity,
using the Z partner's reading frame.

usage: zw_ds.py DATA_DIR OUT_DIR
"""
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.codonalign.codonseq import CodonSeq, _count_diff_NG86, _count_site_NG86, _get_codon_list, cal_dn_ds

TABLE = CodonTable.generic_by_id[1]

data, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(exist_ok=True)
pre = "schistosoma_mansoni.PRJEA36577.WBPS19"

# Strata from smansoni_genome 03_code/smansoni_genome.z_strata.md
AFRICAN = (10_680_112, 23_250_000)
ANCESTRAL = (23_250_001, 43_743_322)


def stratum(chrom, pos):
    if chrom != "SM_V10_Z":
        return "autosome"
    if AFRICAN[0] <= pos <= AFRICAN[1]:
        return "African"
    if ANCESTRAL[0] <= pos <= ANCESTRAL[1]:
        return "ancestral"
    return "PAR"


# --- annotation -------------------------------------------------------------
genes = {}  # gene -> (chrom, start, end, strand, biotype, description)
exons = defaultdict(list)  # transcript -> [(chrom, start, end, strand)]
with gzip.open(data / f"{pre}.annotations.gff3.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            gid = attrs["ID"].removeprefix("gene:")
            genes[gid] = (f[0], int(f[3]), int(f[4]), f[6],
                          attrs.get("biotype", ""), attrs.get("description", "").replace("%2C", ","))
        elif f[2] == "exon":
            for p in attrs["Parent"].split(","):
                exons[p.removeprefix("transcript:")].append((f[0], int(f[3]), int(f[4]), f[6]))


def longest(path):
    best = {}
    with gzip.open(path, "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            g = r.id.split(".")[0]
            if g not in best or len(r.seq) > len(best[g].seq):
                best[g] = r
    return best


prot = longest(data / f"{pre}.protein.fa.gz")
# CDS of the same transcript as the chosen protein
with gzip.open(data / f"{pre}.CDS_transcripts.fa.gz", "rt") as fh:
    _cds = {r.id: r for r in SeqIO.parse(fh, "fasta")}
cds = {g: _cds[r.id] for g, r in prot.items() if r.id in _cds}

w_genes = sorted(g for g, v in genes.items() if v[0] == "SM_V10_WSR")
w_coding = [g for g in w_genes if genes[g][4] == "protein_coding"]
w_pseudo = [g for g in w_genes if genes[g][4] != "protein_coding"]
others = [g for g in prot if genes.get(g, ("",))[0] != "SM_V10_WSR"]


# --- candidate search: shared k-mers, then local alignment ----------------------
def kmers(s, k):
    s = str(s).upper()
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def top_candidates(query_kmers, index, n=15):
    hits = Counter()
    for km in query_kmers:
        for g in index.get(km, ()):
            hits[g] += 1
    return [g for g, _ in hits.most_common(n)]


def build_index(seqs, k):
    idx = defaultdict(set)
    for g, s in seqs.items():
        for km in kmers(s, k):
            idx[km].add(g)
    return idx


blosum = substitution_matrices.load("BLOSUM62")
local_aa = PairwiseAligner(mode="local", substitution_matrix=blosum,
                           open_gap_score=-11, extend_gap_score=-1)
global_aa = PairwiseAligner(mode="global", substitution_matrix=blosum,
                            open_gap_score=-11, extend_gap_score=-1,
                            end_gap_score=0)
global_nt = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                            open_gap_score=-6, extend_gap_score=-2, end_gap_score=0)
local_nt = PairwiseAligner(mode="local", match_score=2, mismatch_score=-3,
                           open_gap_score=-6, extend_gap_score=-2)


def aligned_pairs(aln):
    """Yield (i, j) index pairs of aligned (non-gap) positions."""
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            yield a0 + k, b0 + k


def codon_strings(w_nt, z_nt, codon_map):
    """Build gap-free codon strings from {z_codon_index: w_codon_string}."""
    wc, zc = [], []
    for zi in sorted(codon_map):
        w = codon_map[zi]
        z = str(z_nt[3 * zi:3 * zi + 3]).upper()
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
        except Exception as e:  # saturation / log of negative
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
    diffs = sum(x != y for x, y in zip(w_str, z_str))
    res["nt_identity"] = 1 - diffs / len(w_str)
    return res


rows = []
# W genes whose N-terminus is known to be divergent: also report dN/dS
# either side of this W amino-acid position (U2AF65; Fig. 6 of the paper).
SPLIT = {"Smp_348790": 125}

# --- protein-coding W genes ----------------------------------------------------
aa_idx = build_index({g: prot[g].seq for g in others}, 4)
for w in w_coding:
    wp = str(prot[w].seq).rstrip("*")
    cands = top_candidates(kmers(wp, 4), aa_idx)
    scored = sorted(((local_aa.score(wp, str(prot[z].seq).rstrip("*")), z) for z in cands), reverse=True)
    # Multi-copy families (e.g. histones) give several near-equal protein
    # hits; among those, the gametologue should be the closest in nucleotides.
    near = [z for sc, z in scored if sc >= 0.95 * scored[0][0]]
    w_nt = str(cds[w].seq).upper()
    by_nt = []
    for z in near:
        zp = str(prot[z].seq).rstrip("*")
        pairs = list(aligned_pairs(global_aa.align(wp, zp)[0]))
        z_nt = str(cds[z].seq).upper()
        w_str, z_str = codon_strings(w_nt, z_nt, {j: w_nt[3 * i:3 * i + 3] for i, j in pairs})
        nt_id = sum(a == b for a, b in zip(w_str, z_str)) / max(len(w_str), 1)
        by_nt.append((nt_id, z, pairs, w_str, z_str))
    by_nt.sort(key=lambda x: x[0], reverse=True)
    _, z, pairs, w_str, z_str = by_nt[0]
    zp = str(prot[z].seq).rstrip("*")
    ident = sum(wp[i] == zp[j] for i, j in pairs) / max(len(pairs), 1)
    others_scored = [(sc, g) for sc, g in scored if g != z]
    second = others_scored[0] if others_scored else (0, "")
    c, s, e = genes[z][:3]
    base = {
        "W_gene": w, "W_type": "protein_coding", "W_start": genes[w][1],
        "Z_gene": z, "Z_chrom": c, "Z_start": s, "stratum": stratum(c, s),
        "Z_description": genes[z][5], "aa_identity": ident,
        "W_aa_coverage": len(pairs) / len(wp), "Z_aa_coverage": len(pairs) / len(zp),
        "score": dict((g, sc) for sc, g in scored)[z], "near_ties": len(near) - 1,
        "next_best": second[1], "next_best_score": second[0],
    }
    rows.append({**base, "region": "full", **estimate(w_str, z_str)})
    if w in SPLIT:
        cut = SPLIT[w]
        for region, keep in ((f"aa1-{cut}", lambda i: i < cut), (f"aa{cut + 1}-end", lambda i: i >= cut)):
            sub = {j: w_nt[3 * i:3 * i + 3] for i, j in pairs if keep(i)}
            rows.append({**base, "region": region,
                         **estimate(*codon_strings(w_nt, str(cds[z].seq).upper(), sub))})

# --- W pseudogenes: align spliced exons to candidate CDSs in nucleotide space ---
genome = None
gpath = data / f"{pre}.genomic.fa.gz"
if w_pseudo and gpath.exists():
    try:
        with gzip.open(gpath, "rt") as fh:
            genome = {r.id: r.seq for r in SeqIO.parse(fh, "fasta") if r.id == "SM_V10_WSR"}
    except (EOFError, OSError):  # partial download
        genome = None
if genome:
    nt_idx = build_index({g: cds[g].seq for g in others if g in cds}, 12)
    for w in w_pseudo:
        tx = [t for t in exons if t.split(".")[0] == w]
        ex = sorted(exons[tx[0]], key=lambda x: x[1])
        seq = Seq("".join(str(genome["SM_V10_WSR"][a - 1:b]) for _, a, b, _ in ex))
        if ex[0][3] == "-":
            seq = seq.reverse_complement()
        w_nt = str(seq).upper()
        cands = top_candidates(kmers(w_nt, 12), nt_idx)
        scored = sorted(((local_nt.score(w_nt, str(cds[z].seq).upper()), z) for z in cands), reverse=True)
        best_score, z = scored[0]
        z_nt = str(cds[z].seq).upper()
        aln = local_nt.align(w_nt, z_nt)[0]
        # keep only Z codons whose three bases are all aligned, ungapped
        zpos = {j: i for i, j in aligned_pairs(aln)}
        cmap = {}
        for zi in range(len(z_nt) // 3):
            idx = [zpos.get(3 * zi + k) for k in range(3)]
            if None not in idx and idx[2] - idx[0] == 2:
                cmap[zi] = w_nt[idx[0]:idx[0] + 3]
        w_str, z_str = codon_strings(w_nt, z_nt, cmap)
        c, s, e = genes[z][:3]
        rows.append({
            "W_gene": w, "W_type": "pseudogene", "W_start": genes[w][1],
            "Z_gene": z, "Z_chrom": c, "Z_start": s, "stratum": stratum(c, s),
            "Z_description": genes[z][5], "aa_identity": float("nan"),
            "W_aa_coverage": float("nan"), "Z_aa_coverage": len(cmap) / (len(z_nt) // 3), "region": "full",
            "score": best_score, "next_best": scored[1][1] if len(scored) > 1 else "",
            "next_best_score": scored[1][0] if len(scored) > 1 else 0,
            **estimate(w_str, z_str),
        })
else:
    print(f"skipping {len(w_pseudo)} W pseudogenes: no genome at {gpath}", file=sys.stderr)

df = pd.DataFrame(rows).sort_values("Z_start")
df["dNdS_YN00"] = df["dN_YN00"] / df["dS_YN00"]
df["dNdS_ML"] = df["dN_ML"] / df["dS_ML"]
df.to_csv(out / "zw_gametologue_ds.tsv", sep="\t", index=False, float_format="%.4g")
print(df[["W_gene", "W_type", "region", "Z_gene", "Z_start", "stratum", "aa_identity", "near_ties",
          "codons", "dS_NG86", "dS_YN00", "dS_ML", "dN_ML", "dNdS_ML"]].to_string(index=False))
