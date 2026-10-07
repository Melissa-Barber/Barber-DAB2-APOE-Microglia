#!/usr/bin/env python3
"""
Original AD-Atlas CERAD (high vs low) RELN-core enrichment heatmap,
regenerated from the authoritative fgseaMultilevel results
(ADatlas_CERAD_RELNcore_fgsea_multilevel.csv, "final_CORRECT").

Rows   : EC, ITG, PFC, V1
Columns: RELN core + DAB1 arm, RELN core + DAB2 arm
Cell   : NES (top), significance stars (mid), RELN/APOE-axis leading-edge genes (bottom)
Colour : diverging blue-white-red on NES (-2.5 .. +2.5)
Axis-gene colour: orange = APOE-axis, green = RELN-canonical, black = both.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox

CM = 1 / 2.54
HERE = os.path.dirname(os.path.abspath(__file__))
CSV_CANDS = [
    os.path.join(HERE, "ADatlas_CERAD_RELNcore_fgsea_multilevel.csv"),
    "/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal/Wald-stats multilevel fgsea final_CORRECT/ADatlas_CERAD_RELNcore_fgsea_multilevel.csv",
]
CSV = next((p for p in CSV_CANDS if os.path.exists(p)), CSV_CANDS[0])
df = pd.read_csv(CSV)

REGIONS   = ["EC", "ITG", "PFC", "V1"]
REGION_ID = {r: f"MG_{r}_CERAD_high_vs_low" for r in REGIONS}
ARMS = [("RELN_core_DAB1_v2", "RELN core\n+ DAB1 arm"),
        ("RELN_core_DAB2_v2", "RELN core\n+ DAB2 arm")]

DISPLAY_ORDER = ["APOE", "ITGB1", "RELN"]
GENE_COLOR = {"APOE": "#F2A900",   # APOE-axis  -> orange
              "RELN": "#2CA02C",   # RELN-canonical -> green
              "ITGB1": "#111111"}  # both -> black

cmap = LinearSegmentedColormap.from_list("bwr2", ["#2166AC", "#FFFFFF", "#D6604D"])
norm = Normalize(-2.5, 2.5)
def stars(p): return "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))

STROKE = [pe.withStroke(linewidth=1.8, foreground="white")]
def gene_label(ax, x, y, genes):
    kids = []
    for i, g in enumerate(genes):
        kids.append(TextArea(g, textprops=dict(color=GENE_COLOR.get(g, "#111"),
                     fontsize=8.5, fontweight="bold", path_effects=STROKE)))
        if i < len(genes) - 1:
            kids.append(TextArea(",", textprops=dict(color="#333", fontsize=8.5,
                        fontweight="bold", path_effects=STROKE)))
    pack = HPacker(children=kids, align="center", pad=0, sep=2)
    ax.add_artist(AnnotationBbox(pack, (x, y), frameon=False,
                  box_alignment=(0.5, 0.5), xycoords="data"))

plt.rcParams.update({"font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(18*CM, 15*CM))
fig.suptitle("AD Atlas — CERAD (High vs Low)\nRELN-core enrichment + RELN/APOE-axis leading-edge genes",
             fontsize=11.5, weight="bold", y=1.03)

nrow, ncol = len(REGIONS), len(ARMS)
for ci, (arm_key, _) in enumerate(ARMS):
    for ri, region in enumerate(REGIONS):
        row = df[(df.contrast_id == REGION_ID[region]) & (df.pathway == arm_key)]
        if row.empty:
            continue
        r = row.iloc[0]
        nes = float(r.NES)
        padj = float(r["padj_BH_perContrast"]) if "padj_BH_perContrast" in r else float(r["padj"])
        le = set(str(r.leadingEdge).split(";"))
        genes = [g for g in DISPLAY_ORDER if g in le]
        y = nrow - 1 - ri
        ax.add_patch(plt.Rectangle((ci, y), 1, 1, facecolor=cmap(norm(nes)),
                                   edgecolor="white", lw=3))
        tcol = "white" if nes > 1.1 else "#222"
        ax.text(ci+0.5, y+0.72, f"{nes:+.2f}", ha="center", va="center",
                fontsize=14, weight="bold", color=tcol)
        ax.text(ci+0.5, y+0.52, stars(padj), ha="center", va="center",
                fontsize=12, color=tcol)
        gene_label(ax, ci+0.5, y+0.24, genes)

ax.set_xlim(0, ncol); ax.set_ylim(0, nrow)  # no equal aspect -> wider cells for 3-gene labels
ax.set_xticks([c+0.5 for c in range(ncol)])
ax.set_xticklabels([lab for _, lab in ARMS], fontsize=10.5, weight="bold")
ax.set_yticks([nrow-0.5-i for i in range(nrow)])
ax.set_yticklabels(REGIONS, fontsize=12.5, weight="bold")
ax.xaxis.tick_top(); ax.tick_params(length=0, pad=6)
for s in ax.spines.values():
    s.set_visible(False)

sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03, ticks=[-2, -1, 0, 1, 2])
cbar.ax.set_yticklabels(["≤ -2", "-1", "0", "+1", "≥ +2"], fontsize=8)
cbar.set_label("NES", fontsize=10)

fig.text(0.5, 0.02,
         "v2 sets + Wald-stat (fgseaMultilevel). NES (top), stars *q<0.05 **q<0.01 ***q<0.001; RELN/APOE-axis genes in the leading edge.\n"
         "Axis-gene colour: orange = APOE-axis, green = RELN-canonical, black = both.",
         ha="center", va="bottom", fontsize=7.2, color="#333", style="italic")

fig.subplots_adjust(top=0.82, bottom=0.10, left=0.14, right=0.86)
out = os.path.join(HERE, "Fig_ADatlas_CERAD_RELNcore_leadingEdge.png")
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", out)
