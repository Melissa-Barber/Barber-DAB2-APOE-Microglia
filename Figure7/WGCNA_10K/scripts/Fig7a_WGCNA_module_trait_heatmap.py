#!/usr/bin/env python3
"""
WGCNA top-10,000-gene run (n=39 BIONMG iPSC-microglia, APOE2-excluded) —
two-panel composite re-creating the 5k 'module-trait + gene->module' layout.

  Panel A  Module-trait correlation heatmap
      rows = WGCNA modules ; cols = traits in 3 blocks
        Genotype  : ε4 vs ε3 (Ctrl-siRNA only) | ε4 vs ε3 (all)
        DAB1 KD   : DAB1 (pooled) | 893 | 894
        DAB2 KD   : DAB2 (pooled) | 96  | 98
      cell = point-biserial / Pearson r (blue-white-red, ±1) + BH-FDR stars.

  Panel B  RELN/APOE-axis gene -> module membership bars
      horizontal bars = signed module membership (kME, gene-eigengene r),
      coloured by the gene's module; diamond marker + bold = axis gene.

MAIN : significant module rows only (>=1 BH-sig trait; grey excluded);
       Panel B = axis genes that are members of those significant modules
       (LRP8, LRP1 -> blue ; ITGB8 -> brown ; DAB2 -> purple).
SUPP : all modules; Panel B = all 9 RELN/APOE-axis genes (incl. grey/unassigned
       DAB1 & APOE).

Data (Final Figures /WGCNA_10K/tables):
  module_trait_correlations_expanded_10K_p12.csv  (module x trait: r,padj,sig,n)
  WGCNA_10K_master_gene_table_p12.csv             (kME_own_module, RELN_APOE_axis)
  WGCNA_10K_module_summary_p12.csv                (n_genes, top GO/KEGG)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import matplotlib.patheffects as pe
from scipy import stats as sps

# ---------------------------------------------------------------- style ----
CM = 1 / 2.54
DPI = 600
NES_CMAP = LinearSegmentedColormap.from_list(
    "div_bwr", ["#2166AC", "#FFFFFF", "#B2182B"], N=256)
FONTS = {"title": 14, "row": 13, "col": 12, "value": 9, "star": 12,
         "cbar_lbl": 12, "cbar_tick": 10, "caption": 9.5}
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10, "axes.linewidth": 0.6,
    "savefig.dpi": DPI, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42})
F = FONTS

# ---------------------------------------------------------------- paths ----
ROOT = Path("<EDIT_PROJECT_ROOT>/Final Figures /WGCNA_10K")
TB = ROOT / "tables"
OUT_DIR = ROOT / "figures"; OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = ROOT / "figure_data"; DATA_DIR.mkdir(parents=True, exist_ok=True)
MT_CSV = TB / "module_trait_correlations_expanded_10K_p12.csv"
GENE_CSV = TB / "WGCNA_10K_master_gene_table_p12.csv"
SUM_CSV = TB / "WGCNA_10K_module_summary_p12.csv"
N_SAMPLES = 39

# ---------------------------------------------------------------- colours --
MOD_COL = {
    "turquoise": "#1FB6B6", "blue": "#2C6FB3", "brown": "#7F4F2D",
    "yellow": "#C9A800", "green": "#2E9B3C", "red": "#D62728",
    "black": "#333333", "pink": "#E377A8", "magenta": "#C026B3",
    "purple": "#7E3FAE", "greenyellow": "#7FA300", "tan": "#B08A4F",
    "salmon": "#E8735E", "cyan": "#159AA8", "grey": "#8A8A8A",
}
# short two-line biology annotation (from module summary top GO/KEGG + hubs)
MOD_ANNO = {
    "brown": "Metabolic /\nviral carcinog.",
    "green": "Circadian /\ngene expression",
    "magenta": "Immune resp.\n(FcγR·complement)",
    "blue": "Ribosome /\nchromosome org.",
    "greenyellow": "Lysosome / lipid\n(PPT1·ASAH1)",
    "salmon": "Purinergic /\nER processing",
    "purple": "Lipid metab.\n(DAB2·TREM2)",
    "turquoise": "Tube morphog. /\ncytoskeleton",
    "yellow": "Nucleosome /\nchromatin",
    "pink": "Inflammatory\n(NOD-like)",
    "red": "Cation transport",
    "tan": "Immune resp.\n(TLR8)",
    "cyan": "(MERTK·LRP5)",
    "black": "Antiviral / IFN\n(TLR signalling)",
    "grey": "Unassigned",
}

# Canonical microglial-state per module (ORA marker-overlap; 10K run).
MOD_STATE = {
    "brown": "CRM", "green": "no state", "magenta": "HM", "blue": "RM/HLA",
    "greenyellow": "DAM", "salmon": "HM", "purple": "DAM/LDAM",
}

# ---------------------------------------------------------------- traits ---
# (key, group, label)
TRAITS = [
    ("E4_Ctrl",   "Genotype", "ε4 vs ε3\n(Ctrl)"),
    ("DAB1siRNA", "DAB1 KD",  "DAB1\n(pool)"),
    ("siRNA_893", "DAB1 KD",  "893"),
    ("siRNA_894", "DAB1 KD",  "894"),
    ("DAB2siRNA", "DAB2 KD",  "DAB2\n(pool)"),
    ("siRNA_96",  "DAB2 KD",  "96"),
    ("siRNA_98",  "DAB2 KD",  "98"),
]
TKEY = [t[0] for t in TRAITS]
AXIS_GENES = ["RELN", "LRP8", "LRP1", "ITGB8", "DAB1", "APOE", "DAB2",
              "ITGB1", "VLDLR"]


def stars_padj(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def mm_p(mm, n=N_SAMPLES):
    if pd.isna(mm):
        return np.nan
    r = max(min(float(mm), 0.999999), -0.999999)
    t = r * np.sqrt((n - 2) / (1 - r * r))
    return 2 * sps.t.sf(abs(t), df=n - 2)


# ---------------------------------------------------------------- load -----
mt = pd.read_csv(MT_CSV)
R_tab = mt.pivot(index="module", columns="trait", values="r")
P_tab = mt.pivot(index="module", columns="trait", values="padj")
SIG_tab = mt.pivot(index="module", columns="trait", values="sig").fillna("")

gene = pd.read_csv(GENE_CSV)
summ = pd.read_csv(SUM_CSV).set_index("module")
MOD_SIZE = gene["module"].value_counts().to_dict()  # includes grey

# axis-gene membership
gene_rows = {}
for g in AXIS_GENES:
    s = gene[gene.gene_symbol == g]
    if s.empty:
        gene_rows[g] = dict(gene=g, module=None, MM=np.nan, p=np.nan)
        continue
    r = s.iloc[0]
    gene_rows[g] = dict(gene=g, module=str(r["module"]),
                        MM=float(r["kME_own_module"]),
                        p=mm_p(float(r["kME_own_module"])))

# significant non-grey modules (>=1 BH-sig trait)
sig_mod = sorted({m for m in R_tab.index
                  if m != "grey"
                  and any(SIG_tab.loc[m, t] in ("*", "**", "***") for t in TKEY)},
                 key=lambda m: -R_tab.loc[m, "DAB2siRNA"])
print("trait-significant modules (DAB2-r desc):", sig_mod)

# export the matrix actually plotted
mt.to_csv(DATA_DIR / "Fig1_modtrait_matrix_10K.csv", index=False)
pd.DataFrame(gene_rows.values()).to_csv(
    DATA_DIR / "Fig1_axisgene_membership_10K.csv", index=False)

ALL_MODS = sorted(R_tab.index.tolist(),
                  key=lambda m: (m == "grey", -R_tab.loc[m, "DAB2siRNA"]))

NORM = TwoSlopeNorm(vcenter=0.0, vmin=-1.0, vmax=1.0)


def render(mode):
    keep = sig_mod if mode == "main" else ALL_MODS
    nm = len(keep)
    # MAIN: show only trait columns that have >=1 BH-significant correlation
    # among the displayed (significant) modules; SUPP keeps the full matrix.
    if mode == "main":
        traits = [t for t in TRAITS
                  if any(SIG_tab.loc[m, t[0]] in ("*", "**", "***")
                         for m in keep)]
    else:
        traits = list(TRAITS)
    n_tr = len(traits)
    print(f"[{mode}] trait columns kept: {[t[2].replace(chr(10),' ') for t in traits]}")

    # Panel B gene list
    if mode == "main":
        bars = [(g, gene_rows[g]["module"] or "grey", gene_rows[g]["MM"])
                for g in AXIS_GENES
                if gene_rows[g]["module"] in sig_mod or g == "APOE"]
    else:
        bars = [(g, gene_rows[g]["module"] or "grey",
                 gene_rows[g]["MM"]) for g in AXIS_GENES]
    bars = sorted(bars, key=lambda x: (-(x[2] if not pd.isna(x[2]) else -9)))
    n_bars = len(bars)

    CELL_W, CELL_H = 1.55, 1.55
    XCt = [j * CELL_W for j in range(n_tr)]
    YCm = [i * CELL_H for i in range(nm)]

    BAR_H_UNIT = 0.78
    panelB_h_cm = max(nm * CELL_H, n_bars * BAR_H_UNIT)
    panelA_h_cm = nm * CELL_H

    fig_w_cm = 7.0 + n_tr * CELL_W + 11.5
    fig_h_cm = 6.8 + max(panelA_h_cm, panelB_h_cm)
    fig = plt.figure(figsize=(fig_w_cm * CM, fig_h_cm * CM))

    top = 1 - 2.7 / fig_h_cm
    bottom = 3.6 / fig_h_cm
    band = top - bottom
    left = 0.150
    a_w = (n_tr * CELL_W) / fig_w_cm
    cb_x = left + a_w + 0.012
    b_x = cb_x + 0.11
    b_w = 0.975 - b_x

    aH = band * (panelA_h_cm / max(panelA_h_cm, panelB_h_cm))
    bH = band * (panelB_h_cm / max(panelA_h_cm, panelB_h_cm))
    axA = fig.add_axes([left, top - aH, a_w, aH])
    axB = fig.add_axes([b_x, top - bH, b_w, bH])

    # ---------------- Panel A ----------------
    for i, mod in enumerate(keep):
        for j, t in enumerate(traits):
            r = R_tab.loc[mod, t[0]]; p = P_tab.loc[mod, t[0]]
            xc, yc = XCt[j], YCm[i]
            if pd.isna(r):
                axA.add_patch(Rectangle((xc - CELL_W/2, yc - CELL_H/2),
                              CELL_W, CELL_H, facecolor="white",
                              edgecolor="#CFCFCF", lw=0.6)); continue
            axA.add_patch(Rectangle((xc - CELL_W/2, yc - CELL_H/2),
                          CELL_W, CELL_H, facecolor=NES_CMAP(NORM(r)),
                          edgecolor="white", lw=1.0))
            s = stars_padj(p)
            tcol = "white" if abs(r) > 0.6 else "#111"
            if s:
                axA.text(xc, yc + 0.30, f"{r:+.2f}", ha="center", va="center",
                         fontsize=F["value"], color=tcol, fontweight="bold")
                axA.text(xc, yc - 0.34, s, ha="center", va="center",
                         fontsize=F["star"], color=tcol, fontweight="bold")
                axA.add_patch(Rectangle((xc - CELL_W/2, yc - CELL_H/2),
                              CELL_W, CELL_H, facecolor="none",
                              edgecolor="black", lw=1.4, zorder=6))
            else:
                axA.text(xc, yc, f"{r:+.2f}", ha="center", va="center",
                         fontsize=F["value"], color=tcol)
    axA.set_xlim(-CELL_W/2 - 0.1, (n_tr-1)*CELL_W + CELL_W/2 + 0.1)
    axA.set_ylim((nm-1)*CELL_H + CELL_H/2 + 0.1, -CELL_H/2 - 0.1)
    axA.set_xticks(XCt)
    axA.set_xticklabels([t[2] for t in traits], fontsize=F["col"]-3.5,
                        linespacing=0.95)
    axA.tick_params(axis="x", length=0, pad=4)
    axA.set_yticks(YCm); axA.set_yticklabels([])
    axA.tick_params(axis="y", length=0)
    for sp in axA.spines.values():
        sp.set_visible(False)
    for i, mod in enumerate(keep):
        axA.text(-CELL_W/2 - 0.25, YCm[i],
                 f"{mod} — {MOD_STATE.get(mod,'')} (n={MOD_SIZE[mod]})\n{MOD_ANNO.get(mod,'')}",
                 ha="right", va="center", fontsize=F["row"]-3.0,
                 fontweight="bold", color=MOD_COL[mod], linespacing=0.98)
    # column group brackets
    groups = []
    for j, t in enumerate(traits):
        if not groups or groups[-1][0] != t[1]:
            groups.append([t[1], j, j])
        else:
            groups[-1][2] = j
    gtr = axA.get_xaxis_transform()
    for name, j0, j1 in groups:
        x0 = XCt[j0] - CELL_W/2 + 0.06
        x1 = XCt[j1] + CELL_W/2 - 0.06
        axA.plot([x0, x1], [-0.175, -0.175], transform=gtr, color="#444",
                 lw=1.1, clip_on=False)
        axA.text((XCt[j0] + XCt[j1]) / 2, -0.225, name, transform=gtr,
                 ha="center", va="top", fontsize=F["col"]-1.0,
                 fontweight="bold", color="#333", clip_on=False)

    cax = fig.add_axes([cb_x, top - aH + aH*0.12, 0.012, aH*0.72])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=NES_CMAP, norm=NORM,
                                   ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label("Pearson / point-biserial r", fontsize=F["cbar_lbl"]-2,
                 fontweight="bold")
    cb.ax.tick_params(labelsize=F["cbar_tick"]-1); cb.outline.set_linewidth(0.4)

    # ---------------- Panel B ----------------
    TOPPAD = 1.4
    axB.set_xlim(0.0, 1.12)
    axB.set_ylim(n_bars - 0.5, -0.5 - TOPPAD)
    val_fx = [pe.withStroke(linewidth=1.6, foreground="white")]
    for k, (g, mod, mm) in enumerate(bars):
        col = MOD_COL.get(mod, "#9E9E9E")
        if pd.isna(mm):
            axB.text(0.02, k, f"{g}  —  not in network", ha="left",
                     va="center", fontsize=F["value"], color="#777",
                     style="italic"); continue
        is_grey = (mod == "grey")
        axB.barh(k, mm, height=0.62, color=col, edgecolor="black",
                 linewidth=1.4, zorder=3, alpha=0.55 if is_grey else 1.0)
        if abs(mm) > 0.16:
            axB.text(mm / 2, k, f"{mm:+.2f}", ha="center", va="center",
                     fontsize=F["value"]-0.5, color="#111", zorder=6,
                     path_effects=val_fx)
        mk_x, tx_x = mm + 0.012, mm + 0.030
        axB.plot(mk_x, k, marker="D", ms=5.0, color="#111", zorder=5,
                 clip_on=False)
        lab = g + ("  (grey/unassigned)" if is_grey else "")
        axB.text(tx_x, k, lab, ha="left", va="center", fontsize=F["value"]+1.5,
                 fontweight="bold", color="#111", zorder=4)
    axB.axvline(0, color="#888", lw=0.8)
    for xthr in (0.6, 0.8):
        axB.axvline(xthr, color="#555", lw=1.0, linestyle="--", zorder=2)
    for xmid, lab in [(0.70, "peripheral\n0.6–0.8"), (0.93, "hub\n≥0.8")]:
        axB.text(xmid, -0.5 - TOPPAD * 0.55, lab, ha="center", va="center",
                 fontsize=F["value"]-1.0, color="#555", fontweight="bold",
                 linespacing=0.95)
    axB.set_yticks([])
    for sp in ("top", "right", "left"):
        axB.spines[sp].set_visible(False)
    axB.set_xlabel("Module membership (signed kME = gene–eigengene r)",
                   fontsize=F["col"]-1.5, labelpad=7)
    axB.tick_params(axis="x", labelsize=F["col"]-2)
    axB.grid(axis="x", color="#e3e3e3", lw=0.5, zorder=0)

    # ---- headers ----
    def header(ax, text):
        pos = ax.get_position()
        bh = 0.040; by = pos.y1 + 0.014
        fig.patches.append(FancyBboxPatch((pos.x0, by), pos.width, bh,
                           boxstyle="round,pad=0.004,rounding_size=0.012",
                           transform=fig.transFigure, facecolor="#EBEBEB",
                           edgecolor="none", clip_on=False, zorder=10))
        fig.text(pos.x0 + pos.width/2, by + bh/2, text, ha="center",
                 va="center", fontsize=F["title"]-1, fontweight="bold",
                 color="#222", zorder=11)
    header(axA, "A.  Module–trait correlation")
    header(axB, "B.  RELN/APOE-axis members"
                if mode == "main" else "B.  RELN/APOE-axis genes → module")

    ttl = ("WGCNA (top-10,000-gene run) module–trait correlation & "
           "RELN/APOE-axis membership — n=39 BIONMG iPSC-MG")
    fig.suptitle(ttl + ("" if mode == "main" else " — supplementary"),
                 x=0.012, y=0.985, ha="left",
                 fontsize=F["title"], fontweight="bold")

    foot = (
        "Panel A: point-biserial / Pearson r; black border + stars = BH-FDR sig "
        "(***padj<0.001, **<0.01, *<0.05).  ε4 vs ε3 (Ctrl)=genotype in "
        "control-siRNA only (n=8).  DAB1 KD: pooled (n=25), "
        "893 (n=17), 894 (n=16) vs ctrl.  DAB2 KD: pooled (n=22), 96 (n=14), "
        "98 (n=16) vs ctrl.  Panel B: diamond + bold = RELN/APOE-axis gene; "
        "kME computed across all 39 samples."
    )
    if mode == "main":
        foot += ("  MAIN: significant module rows AND significant trait columns "
                 "only — a trait is shown only if ≥1 displayed module reaches "
                 "BH-padj<0.05 (ε4-vs-ε3, DAB1-pool and 894 had none and are "
                 "dropped; see supplementary for the full matrix).  Panel B = "
                 "axis members of those modules (LRP8, LRP1→blue; ITGB8→brown; "
                 "DAB2→purple) plus APOE (grey/unassigned, kME≈0).")
    else:
        foot += ("  SUPP: all 14 modules + grey; Panel B = all 9 RELN/APOE-axis "
                 "genes (DAB1 & APOE fall in grey/unassigned).")
    fig.text(0.012, 0.012, foot, ha="left", va="bottom",
             fontsize=F["caption"]-2.0, color="#555", style="italic", wrap=True)

    out = OUT_DIR / (f"Fig1_10K_modtrait_membership_"
                     f"{'main' if mode == 'main' else 'supp'}.png")
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"saved {out.name}  (modules={nm}, barsB={n_bars})")


if __name__ == "__main__":
    render("main")
    render("supp")
