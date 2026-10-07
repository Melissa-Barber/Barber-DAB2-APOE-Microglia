#!/usr/bin/env python3
"""
FigS4c — Hammond 2019 developmental MG clusters × RELN_core gene overlap.

Re-rendered in the manuscript-standard palette: a sequential white →
#B2182B (the warm half of the diverging blue-white-red used elsewhere)
so this percent-overlap heatmap reads as visually consistent with the
NES / Pearson-r heatmaps in the same figure family.

Cell content: overlap count, % of RELN set, significance stars.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import apply_style, FONTS, CM, DPI
apply_style()

CSV     = Path(__file__).resolve().parents[1] / "data" / "FigS4c_RELN_Hammond_overlap.csv"
OUT_DIR = Path(__file__).resolve().parents[1] / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT     = OUT_DIR / "FigS4c_Hammond_RELN_overlap_heatmap.png"

# Sequential palette = pale cream → manuscript red
SEQ_CMAP = LinearSegmentedColormap.from_list(
    "white_to_red", ["#FFFFFF", "#FBEAE5", "#F4B49E", "#D77B7B", "#B2182B"], N=256)

ROWS = [
    ("RELN_core_DAB1", "RELN core\n+ DAB1 arm\n(n=115)"),
    ("RELN_core_DAB2", "RELN core\n+ DAB2 arm\n(n=113)"),
]
# Display label includes the cluster's gene-set size in parentheses on a
# second line so the reader can put each overlap count and percentage in
# context. Gene-set sizes are taken from mg_hammond_li2019_clusters_long.csv.
COLS_DISPLAY = [
    ("C1 (Arg1+/early fetal)",       "C1 (Arg1+/early fetal)\n(n=98)"),
    ("C2a (Proliferating)",          "C2a (Proliferating)\n(n=100)"),
    ("C3 (Fabp5+/metabolic)",        "C3 (Fabp5+/metabolic)\n(n=99)"),
    ("C4 (Spp1+/ATM)",               "C4 (Spp1+/ATM)\n(n=98)"),
    ("C5 (Hmox1+/stress)",           "C5 (Hmox1+/stress)\n(n=100)"),
    ("C6 (Ms4a7+/border)",           "C6 (Ms4a7+/border)\n(n=98)"),
    ("C7a (Homeostatic)",            "C7a (Homeostatic)\n(n=44)"),
    ("C7b (Homeostatic+complement)", "C7b (Homeostatic\n+ complement)\n(n=90)"),
]
COLS     = [c[1] for c in COLS_DISPLAY]
COL_KEYS = [c[0] for c in COLS_DISPLAY]   # match CSV strings

df = pd.read_csv(CSV)
pct_mat = np.full((len(ROWS), len(COL_KEYS)), np.nan)
n_mat   = np.zeros_like(pct_mat, dtype=int)
sig_mat = np.full_like(pct_mat, "", dtype=object)
for i, (rk, _) in enumerate(ROWS):
    for j, ck in enumerate(COL_KEYS):
        sub = df[(df["reln_set"] == rk) & (df["hammond_cluster"] == ck)]
        if not sub.empty:
            pct_mat[i, j] = float(sub["pct"].iloc[0])
            n_mat[i, j]   = int(sub["overlap"].iloc[0])
            sig_mat[i, j] = str(sub["sig"].iloc[0])

# Map "." (trend) and "ns" so the rendering keeps the user's existing convention
def disp_sig(s):
    if s in ("***", "**", "*"): return s
    return "ns"   # near-significant (padj 0.05-0.10) now shown as ns (dot removed)

vmax = max(5.0, float(np.nanmax(pct_mat)) * 1.5)
norm = Normalize(vmin=0, vmax=vmax)

# ---- Figure --------------------------------------------------------------
fig_w_cm, fig_h_cm = 26.0, 13.5
fig, ax = plt.subplots(figsize=(fig_w_cm * CM, fig_h_cm * CM))
fig.subplots_adjust(left=0.17, right=0.86, top=0.91, bottom=0.46)

CELL_W = 2.0
CELL_H = 1.6
for i in range(len(ROWS)):
    for j in range(len(COLS)):
        xc = j * CELL_W
        yc = i * CELL_H
        v = pct_mat[i, j]
        face = SEQ_CMAP(norm(v)) if not np.isnan(v) else "white"
        ax.add_patch(Rectangle((xc - CELL_W/2, yc - CELL_H/2), CELL_W, CELL_H,
                                  facecolor=face, edgecolor="#666", linewidth=0.8))
        if np.isnan(v):
            continue
        text_col = "white" if v > vmax * 0.55 else "#222"
        s = disp_sig(sig_mat[i, j])
        is_sig = s in ("***", "**", "*")
        ax.text(xc, yc - 0.34, f"{n_mat[i, j]}  ({v:.1f}%)",
                ha="center", va="center",
                fontsize=FONTS["value"] + 3, color=text_col,
                weight=("bold" if is_sig else "normal"))
        ax.text(xc, yc + 0.20, s, ha="center", va="center",
                fontsize=FONTS["star"] + 3, color=text_col,
                weight=("bold" if is_sig else "normal"))

ax.set_xlim(-CELL_W/2 - 0.15,
            (len(COLS) - 1) * CELL_W + CELL_W/2 + 0.15)
ax.set_ylim((len(ROWS) - 1) * CELL_H + CELL_H/2 + 0.15,
            -CELL_H/2 - 0.15)
ax.set_yticks([i * CELL_H for i in range(len(ROWS))])
ax.set_yticklabels([r[1] for r in ROWS], fontsize=FONTS["row"], weight="bold")
ax.set_xticks([j * CELL_W for j in range(len(COLS))])
ax.set_xticklabels(COLS, fontsize=FONTS["col"] - 1, rotation=45, ha="right")
ax.tick_params(axis="x", pad=20, top=False, bottom=False)
ax.tick_params(axis="y", left=False, right=False)
for sp in ("top", "right", "bottom", "left"):
    ax.spines[sp].set_visible(False)

# Title
fig.suptitle("Hammond 2019 developmental MG clusters — gene overlap with RELN_core",
             fontsize=FONTS["title"], weight="bold", y=0.96)

# Colourbar
cax = fig.add_axes([0.88, 0.35, 0.022, 0.45])
cb = mpl.colorbar.ColorbarBase(cax, cmap=SEQ_CMAP, norm=norm,
                                  ticks=np.linspace(0, vmax, 4).round(1))
cb.set_label("% overlap\n(of RELN set)",
             fontsize=FONTS["cbar_lbl"], weight="bold", labelpad=2)
cb.ax.tick_params(labelsize=FONTS["cbar_tick"])

# Footnote describing the test
fig.text(0.50, 0.04,
         "Hypergeometric test, BH-FDR; *** padj<0.001, ** padj<0.01, "
         "* padj<0.05, ns ≥ 0.05.   Universe N = 20 000 genes.",
         ha="center", va="bottom",
         fontsize=FONTS["caption"] - 1, color="#444", style="italic")

fig.savefig(OUT, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"saved {OUT}")
