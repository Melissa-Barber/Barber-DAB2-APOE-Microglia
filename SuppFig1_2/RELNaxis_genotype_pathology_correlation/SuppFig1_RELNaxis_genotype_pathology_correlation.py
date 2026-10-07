#!/usr/bin/env python3
"""
RELN/APOE-axis genes x AD pathology (Composite Stage / Braak / CERAD) per region,
genotype-stratified correlation table (E3-driven vs e4-driven).

Reconstructed standalone renderer, matched to the original compact/wide layout.
The rotated "Canonical / APOE-canonical RELN" row-group text label on the far
left has been removed per request; the coloured left group strip and gene
labels are kept.

Cell = Spearman r (top), genotype badge on driven cells, significance stars (bottom).
Cell fill / text colour: blue = E3/E3-driven, red = e4-driven, purple = both, grey = ns.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

CM = 1 / 2.54
HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(os.path.dirname(HERE), "data",
                   "RELNaxis_genotype_pathology_correlation.csv")
OUT = os.path.join(os.path.dirname(HERE), "figures",
                   "RELNaxis_genotype_pathology_correlation.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

df = pd.read_csv(CSV, keep_default_na=False)

GENES = ["RELN", "DAB1", "LRP8", "ITGB1", "DAB2", "LRP1", "APOE"]
GENE_GROUP = {"RELN": "Canonical", "DAB1": "Canonical", "LRP8": "Canonical", "ITGB1": "Canonical",
              "DAB2": "APOE", "LRP1": "APOE", "APOE": "APOE"}
GENE_COLOR = {"RELN": "#117733", "DAB1": "#C0392B", "LRP8": "#7D3C98", "ITGB1": "#1E8449",
              "DAB2": "#2471A3", "LRP1": "#E07A5F", "APOE": "#C9A227"}
GROUP_STRIP = {"Canonical": "#16416b", "APOE": "#8c2b2b"}
METRICS = ["Composite Stage", "Braak", "CERAD"]
REGIONS = ["EC", "ITG", "V1", "PFC"]

BLUE_FILL, BLUE_TXT = "#DCEBF7", "#1F5FA5"
RED_FILL,  RED_TXT  = "#F9DEDE", "#C0392B"
PUR_FILL,  PUR_TXT  = "#EADCF2", "#7D3C98"
NS_TXT = "#9AA0A6"

def cell_style(geno):
    if geno == "E3":  return BLUE_FILL, BLUE_TXT
    if geno == "e4":  return RED_FILL,  RED_TXT
    if geno == "both":return PUR_FILL,  PUR_TXT
    return "#FFFFFF", NS_TXT

def badge_txt(geno):
    return {"E3": "E3", "e4": "ε4", "both": "E3/ε4"}.get(geno, "")

# ---- geometry (compact/wide: cells wider than tall) -----------------------
CW, CH = 1.0, 0.56
X0 = 1.9                       # width of the Gene column (strip + gene labels)
ncols = len(METRICS) * len(REGIONS)

# fonts (small, matching the original)
F_GENE_HDR, F_METRIC, F_REGION = 10.5, 10.5, 9.5
F_GENE, F_VAL, F_BADGE, F_SIG, F_NS = 10.5, 11.5, 8.5, 10.0, 9.5
F_LEG, F_SIGKEY = 9.5, 8.5

fig, ax = plt.subplots(figsize=(15 * CM * 2.6, 7.2 * CM * 2.6))

def col_x(mi, ri):
    return X0 + (mi * len(REGIONS) + ri) * CW

hdr_y = len(GENES) * CH
HDR = CH * 1.05                # each header sub-row height
# Gene header block
ax.add_patch(Rectangle((0, hdr_y), X0, HDR*2, facecolor="#f0f0f0", edgecolor="#bbb", lw=0.8))
ax.text(X0/2, hdr_y + HDR, "Gene", ha="center", va="center", fontsize=F_GENE_HDR, weight="bold")
for mi, m in enumerate(METRICS):
    gx0 = col_x(mi, 0); gx1 = col_x(mi, len(REGIONS)-1) + CW
    ax.add_patch(Rectangle((gx0, hdr_y+HDR), gx1-gx0, HDR,
                 facecolor="#d9d9d9", edgecolor="#bbb", lw=0.8))
    ax.text((gx0+gx1)/2, hdr_y + HDR*1.5, m, ha="center", va="center",
            fontsize=F_METRIC, weight="bold")
    for ri, rg in enumerate(REGIONS):
        cx = col_x(mi, ri)
        ax.add_patch(Rectangle((cx, hdr_y), CW, HDR,
                     facecolor="#ececec", edgecolor="#bbb", lw=0.8))
        ax.text(cx+CW/2, hdr_y+HDR*0.5, rg, ha="center", va="center",
                fontsize=F_REGION, weight="bold")

# ---- coloured group strips (no rotated text) ------------------------------
for grp, genes in (("Canonical", GENES[:4]), ("APOE", GENES[4:])):
    ys = [len(GENES)-1-GENES.index(g) for g in genes]
    y_lo = min(ys)*CH; y_hi = (max(ys)+1)*CH
    ax.add_patch(Rectangle((0.05, y_lo), 0.20, y_hi-y_lo,
                 facecolor=GROUP_STRIP[grp], edgecolor="none"))

# ---- gene rows ------------------------------------------------------------
for gi, gene in enumerate(GENES):
    y = (len(GENES)-1-gi) * CH
    ax.add_patch(Rectangle((0.25, y), X0-0.25, CH,
                 facecolor=("#f7fbff" if GENE_GROUP[gene]=="Canonical" else "#fdf6f4"),
                 edgecolor="#ddd", lw=0.6))
    ax.text((0.25+X0)/2, y+CH/2, gene, ha="center", va="center",
            fontsize=F_GENE, weight="bold", color=GENE_COLOR[gene])
    for mi, m in enumerate(METRICS):
        for ri, rg in enumerate(REGIONS):
            cx = col_x(mi, ri)
            sub = df[(df.gene==gene)&(df.metric==m)&(df.region==rg)]
            if sub.empty: continue
            r = float(sub.r.iloc[0]); geno = str(sub.geno.iloc[0]); sig = str(sub.sig.iloc[0])
            fill, vtxt = cell_style(geno)
            ax.add_patch(Rectangle((cx, y), CW, CH, facecolor=fill, edgecolor="#cfcfcf", lw=0.7))
            is_sig = geno in ("E3","e4","both")
            vy = 0.49 if is_sig else 0.55   # drop the number lower in badge cells so it clears the E3/e4 tag
            ax.text(cx+CW/2, y+CH*vy, f"{r:+.2f}", ha="center", va="center",
                    fontsize=F_VAL, weight=("bold" if is_sig else "normal"), color=vtxt)
            if is_sig:
                bcol = BLUE_TXT if geno=="E3" else (RED_TXT if geno=="e4" else PUR_TXT)
                ax.add_patch(FancyBboxPatch((cx+CW*0.28, y+CH*0.70), CW*0.44, CH*0.24,
                             boxstyle="round,pad=0.015,rounding_size=0.04",
                             facecolor=bcol, edgecolor="none", mutation_aspect=1))
                ax.text(cx+CW/2, y+CH*0.82, badge_txt(geno), ha="center", va="center",
                        fontsize=F_BADGE, weight="bold", color="white")
                ax.text(cx+CW/2, y+CH*0.15, sig, ha="center", va="center",
                        fontsize=F_SIG, color=vtxt, weight="bold")
            else:
                ax.text(cx+CW/2, y+CH*0.24, "ns", ha="center", va="center",
                        fontsize=F_NS, color=vtxt, style="italic")

ax.add_patch(Rectangle((0.25, 0), (X0-0.25)+ncols*CW, len(GENES)*CH,
             fill=False, edgecolor="#888", lw=1.2))

# ---- legend + significance key --------------------------------------------
leg_y = -0.95; bh = 0.34
items = [("E3/E3-driven", BLUE_FILL, BLUE_TXT), ("ε4-driven", RED_FILL, RED_TXT),
         ("Both", PUR_FILL, PUR_TXT), ("ns", "#FFFFFF", "#888")]
lx = X0
for lab, fc, ec in items:
    ax.add_patch(FancyBboxPatch((lx, leg_y), bh, bh,
                 boxstyle="round,pad=0.01,rounding_size=0.06",
                 facecolor=fc, edgecolor=ec, lw=1.4))
    ax.text(lx+bh+0.12, leg_y+bh/2, lab, ha="left", va="center", fontsize=F_LEG,
            weight="bold", color=ec)
    lx += bh + 0.20 + 0.115*len(lab)
ax.text(col_x(len(METRICS)-1, len(REGIONS)-1)+CW, leg_y+bh/2,
        "Significance:  *** p<0.001    ** p<0.01    * p<0.05",
        ha="right", va="center", fontsize=F_SIGKEY, color="#333")

ax.set_xlim(-0.05, X0 + ncols*CW + 0.1)
ax.set_ylim(leg_y-0.3, len(GENES)*CH + HDR*2 + 0.08)
ax.set_aspect("equal"); ax.axis("off")
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", OUT)
