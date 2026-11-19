"""Plot Z-W gametologue synonymous divergence along the S. mansoni Z.

usage: plot_ds.py zw_gametologue_ds.tsv OUT_PREFIX
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t")
df = df[(df.region == "full") & (df.stratum != "PAR")].copy()
df["mb"] = df.Z_start / 1e6
CAP = 1.2

fig, axes = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
for ax in axes:
    ax.axvspan(10.680112, 23.25, color="#F1CFED", lw=0)
    ax.axvspan(23.25, 43.743322, color="#FDE6FA", lw=0)
    ax.axvspan(43.743322, 51.3, color="#EEFDE6", lw=0)

panels = [("pS", "pS (proportion synonymous\ndifferences, NG86)", 0.8),
          ("dS_ML", "dS (ML, GY94)", CAP)]
for ax, (col, label, ymax) in zip(axes, panels):
    for kind, marker in (("protein_coding", "o"), ("pseudogene", "^")):
        d = df[df.W_type == kind]
        if d.empty:
            continue
        y = d[col].clip(upper=ymax * 0.97)
        tie = d.near_ties > 0
        ax.scatter(d.mb[~tie], y[~tie], s=18, marker=marker, color="black",
                   label="W " + kind.replace("_", "-"))
        ax.scatter(d.mb[tie], y[tie], s=18, marker=marker, facecolor="none",
                   edgecolor="black", label="ambiguous pair (multi-copy)" if kind == "protein_coding" else None)
    sat = df[df[col] >= ymax * 0.97]
    for _, r in sat.iterrows():
        ax.annotate("saturated", (r.mb, ymax * 0.97), xytext=(4, -3),
                    textcoords="offset points", fontsize=7)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(label, fontsize=8)

u2af = df[df.W_gene == "Smp_348790"].iloc[0]
axes[0].annotate("U2AF65", (u2af.mb, min(u2af.pS, 0.78)), xytext=(4, -12),
                 textcoords="offset points", fontsize=7)
axes[1].legend(fontsize=7, loc="upper left", frameon=False)
axes[1].set_xlabel("Position of Z gametologue on SM_V10_Z (Mb)")
axes[1].set_xlim(10, 45)
axes[0].text(10.9, 0.74, "African stratum", fontsize=7)
axes[0].text(33, 0.74, "ancestral stratum", fontsize=7)
fig.tight_layout()
fig.savefig(sys.argv[2] + ".pdf")
fig.savefig(sys.argv[2] + ".png", dpi=200)
