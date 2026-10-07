#!/usr/bin/env python3
"""
RELN-core + MG states — APOE-ε4 vs ε3 (BIONMG only)

Single-column heatmap: BIONMG iPSC-MG Ctrl-siRNA E4 vs E3
Gene sets: RELN-core-DAB1/2, MG states (Mancuso), Ribosomal,
           LDAM (Haney), APOE4 lipid iMG (Victor)

Run:
    python3 "Figures/Claude revised scripts/APOE_E4vsE3_MGstates_RELNcore_heatmap/scripts/render_heatmap.py"
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))  # style.py bundled alongside this script
from style import apply_style, NES_CMAP, CM, DPI
from matplotlib.colors import TwoSlopeNorm
apply_style()

NES_NORM = TwoSlopeNorm(vcenter=0.0, vmin=-3, vmax=3)

FONTS = {
    "title":     13,
    "subtitle":   9,
    "col":       11,
    "row":       11,
    "value":     13,
    "star":      15,
    "le_info":    9,
    "cbar_lbl":  12,
    "cbar_tick": 10,
}

ROOT    = Path(__file__).resolve().parents[1]
DATA    = ROOT / "data" / "BIONMG_E4vsE3_Ctrl_MGstates_RELNcore_fgsea.csv"
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(DATA)
print(df[["pathway", "NES", "padj", "size"]].to_string())

ROWS = [
    ("RELN-core-DAB1",         "RELN_core_DAB1",              "reln"),
    ("RELN-core-DAB2",         "RELN_core_DAB2",              "reln"),
    ("MG Homeostatic (HM)",    "MG_HM",                       "mg"),
    ("MG DAM",                 "MG_DAM",                      "mg"),
    ("MG HLA",                 "MG_HLA",                      "mg"),
    ("MG CRM",                 "MG_CRM",                      "mg"),
    ("MG IRM",                 "MG_IRM",                      "mg"),
    ("Ribosomal (RM)",         "Ribosomal_RM",                "mg"),
    ("LDAM markers (Haney)",   "LDAM_Haney2024",              "lipid"),
    ("APOE4 lipid iMG (Victor)", "Victor2022_APOE4_Lipid_iMG", "lipid"),
]

COL_LABEL = "BIONMG iMG\nE4 vs E3\n(Ctrl siRNA)"

RELN_AXIS = {"DAB1","DAB2","APOE","LRP8","ITGB1","ITGB8","VLDLR"}


def stars(p):
    if pd.isna(p): return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""


def parse_le(le_str):
    if pd.isna(le_str) or le_str == "":
        return []
    return [g.strip().upper() for g in str(le_str).split(";") if g.strip()]


data_rows = []
for disp, key, block in ROWS:
    sub = df[df["pathway"] == key]
    if sub.empty:
        data_rows.append({
            "label": disp, "key": key, "block": block,
            "NES": np.nan, "stars": "", "le_axis": "", "le_count": 0,
            "size": 0,
        })
        continue
    r = sub.iloc[0]
    le_genes = parse_le(r.get("leadingEdge", ""))
    le_axis = sorted(set(g for g in le_genes if g in RELN_AXIS))
    data_rows.append({
        "label": disp, "key": key, "block": block,
        "NES": float(r["NES"]),
        "stars": stars(float(r["padj"])),
        "le_axis": ", ".join(le_axis) if le_axis else "",
        "le_count": len(le_genes),
        "size": int(r["size"]),
    })

BLOCK_GAP = 0.5
CELL_W = 4.0
CELL_H = 1.6

def y_for(i):
    b = ROWS[i][2]
    if b == "reln":  return i
    if b == "mg":    return i + BLOCK_GAP
    return i + 2 * BLOCK_GAP

n = len(data_rows)
y_max = y_for(n - 1) + 0.6

fig_w_cm = 16.0
fig_h_cm = 6.0 + 1.85 * y_max          # taller rows so gene names are not squashed
fig, ax = plt.subplots(figsize=(fig_w_cm * CM, fig_h_cm * CM))
fig.subplots_adjust(left=0.42, right=0.72, top=0.80, bottom=0.06)

for i, rd in enumerate(data_rows):
    yc = y_for(i) * CELL_H
    xc = 0
    v = rd["NES"]

    if np.isnan(v):
        face = (0.93, 0.93, 0.93, 1.0)
        txt = "—"
    else:
        face = NES_CMAP(NES_NORM(np.clip(v, -3, 3)))
        txt = f"{v:+.2f}"

    ax.add_patch(Rectangle(
        (xc - CELL_W / 2, yc - CELL_H / 2), CELL_W, CELL_H,
        facecolor=face, edgecolor="white", linewidth=1.2,
    ))

    text_col = "white" if (not np.isnan(v) and abs(v) > 1.85) else "#222"

    star_str = rd["stars"]

    le_info = ""
    if rd["le_count"] > 0 and rd["size"] > 0:
        le_info = f"({rd['le_count']}/{rd['size']})"

    value_line = f"{txt}  {le_info}".strip()

    # Stack up to three horizontally-centred zones inside the cell:
    #   1. NES value (+ leading-edge count)      — bold
    #   2. significance star                      — bold, larger, sits mid-cell
    #   3. RELN/APOE-axis leading-edge genes      — italic, its own uncrowded line
    # Giving the star its own centred line keeps it clearly visible at the middle
    # of the cell instead of trailing off the right edge of the gene-name line.
    lines = [("value", value_line)]
    if star_str:
        lines.append(("star", star_str))
    if rd["le_axis"]:
        lines.append(("gene", rd["le_axis"]))

    line_gap = 0.30 * CELL_H
    n_lines = len(lines)
    y0 = yc - (n_lines - 1) / 2.0 * line_gap
    for j, (kind, text) in enumerate(lines):
        yj = y0 + j * line_gap
        if kind == "value":
            ax.text(xc, yj, text, ha="center", va="center",
                    fontsize=FONTS["value"], color=text_col, weight="bold")
        elif kind == "star":
            ax.text(xc, yj, text, ha="center", va="center",
                    fontsize=FONTS["star"], color=text_col, weight="bold")
        else:  # gene
            ax.text(xc, yj, text, ha="center", va="center",
                    fontsize=FONTS["le_info"], color=text_col, style="italic")

# Block separators
for i in range(1, n):
    if ROWS[i][2] != ROWS[i - 1][2]:
        sep_y = (y_for(i) + y_for(i - 1)) / 2 * CELL_H
        ax.plot([-CELL_W / 2 - 0.1, CELL_W / 2 + 0.1],
                [sep_y, sep_y],
                linestyle="--", color="#666", linewidth=0.9)

ax.set_xticks([0])
ax.set_xticklabels([COL_LABEL], fontsize=FONTS["col"], weight="bold")
ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")

ax.set_yticks([y_for(i) * CELL_H for i in range(n)])
ax.set_yticklabels([r["label"] for r in data_rows], fontsize=FONTS["row"])

ax.set_xlim(-CELL_W / 2 - 0.2, CELL_W / 2 + 0.2)
ax.set_ylim(y_max * CELL_H + 0.2, -CELL_H / 2 - 0.4)
ax.tick_params(top=False, bottom=False, left=False, right=False, pad=10)
for sp in ("top", "right", "bottom", "left"):
    ax.spines[sp].set_visible(False)

cax = fig.add_axes([0.78, 0.12, 0.025, 0.54])
cb = mpl.colorbar.ColorbarBase(cax, cmap=NES_CMAP, norm=NES_NORM,
                                ticks=[-3, -2, -1, 0, 1, 2, 3])
cb.set_label("NES", fontsize=FONTS["cbar_lbl"], weight="bold")
cb.ax.tick_params(labelsize=FONTS["cbar_tick"])
cb.outline.set_linewidth(0.4)

fig.text(0.50, 0.965,
         "RELN-core + MG states — APOE-ε4 vs ε3 (BIONMG)",
         ha="center", va="top",
         fontsize=FONTS["title"], weight="bold")
fig.text(0.50, 0.925,
         "R fgseaMultilevel · v2 sets · Wald-stat · "
         "NES shown; * padj<0.05  ** <0.01  *** <0.001",
         ha="center", va="top",
         fontsize=FONTS["subtitle"], color="#555")

out = OUT_DIR / "APOE_E4vsE3_MGstates_RELNcore_BIONMG_only.png"
fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"\nsaved {out}")
