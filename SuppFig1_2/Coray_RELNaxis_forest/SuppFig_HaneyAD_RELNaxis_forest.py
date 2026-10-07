#!/usr/bin/env python3
"""
Coray/Haney iPSC-microglia (APOE4/4) fAβ arm — RELN/APOE-axis gene log2 fold-changes
across two contrasts (AB vs NT; AB+GNE vs AB): DESeq2 shrunken-LFC forest panel.

Renders directly from the exact DESeq2 stats behind the original figure
(Coray_RELNaxis_forest_data.csv), so every value — including DAB2 (ENSG00000153071)
— matches the source. x-axis label: "log2 fold-change (DESeq2)" (", shrunken" removed).
DAB2 role annotation: "RELN and APOE receptors endocytic adaptor protein".
"""
import os, textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

CM, DPI = 1/2.54, 300
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

HERE = os.path.dirname(os.path.abspath(__file__))
CSV  = os.path.join(os.path.dirname(HERE), "data", "Coray_RELNaxis_forest_data.csv")
OUT  = os.path.join(os.path.dirname(HERE), "figures", "Coray_RELNaxis_forest.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

tab = pd.read_csv(CSV)
GENES = list(dict.fromkeys(tab["gene"]))
ROLE  = {g: tab[tab.gene == g]["role"].iloc[0] for g in GENES}
CONTRASTS = [("AB vs NT", "#D7263D"), ("AB + GNE vs AB", "#1F77B4")]

def stars(q):
    if pd.isna(q): return ""
    return "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else ""))

def wrap_role(role, width=26):
    return "\n".join(textwrap.wrap(role, width=width)) if len(role) > width else role

n_genes = len(GENES)
fig, ax = plt.subplots(figsize=(21*CM, 10.8*CM), dpi=DPI)

y_band = np.arange(n_genes)
dy = 0.18
for i_g, sym in enumerate(GENES):
    for j_c, (cname, color) in enumerate(CONTRASTS):
        sub = tab[(tab.gene == sym) & (tab.contrast == cname)]
        if sub.empty: continue
        r = sub.iloc[0]
        lfc, se, q = float(r.log2FC), float(r.lfcSE), float(r.padj)
        if not np.isfinite(lfc): continue
        y = y_band[i_g] + (j_c - 0.5) * dy * 2
        ax.errorbar(lfc, y, xerr=se, fmt="o", color=color, ecolor=color,
                    ms=7, elinewidth=1.4, capsize=4, mec="#222", mew=0.6, zorder=3)
        s = stars(q)
        lbl = f"{lfc:+.2f}{(' '+s) if s else ''}"
        offx = 0.06 if lfc >= 0 else -0.06
        xt = (lfc + se + offx) if lfc >= 0 else (lfc - se + offx)
        ax.text(xt, y, lbl, ha=("left" if lfc >= 0 else "right"),
                va="center", fontsize=8, color="#222")

ax.axvline(0, color="#444", lw=0.8, ls="--", zorder=1)
for i in range(n_genes):
    if i % 2 == 0:
        ax.axhspan(i - 0.5, i + 0.5, color="#F7F7F7", zorder=0)

ax.set_yticks(y_band)
ax.set_yticklabels([f"{sym}\n({wrap_role(ROLE[sym])})" for sym in GENES], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("log₂ fold-change (DESeq2)", fontsize=10)
ax.tick_params(axis="x", labelsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_xlim(-3.0, 3.0)

for cname, color in CONTRASTS:
    ax.plot([], [], "o", color=color, label=cname, ms=7, mec="#222", mew=0.6)
ax.legend(loc="upper right", fontsize=8, frameon=False,
          title="Contrast", title_fontsize=8)

plt.tight_layout()
fig.savefig(OUT, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"saved {OUT}")
