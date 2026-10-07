#!/usr/bin/env python3
"""
WGCNA top-10,000-gene run (n=39 BIONMG iPSC-MG) — module × microglial-state
composite (re-creating the 5k composite_MGstate_moduletrait layout).

Two stacked panels sharing the SAME module columns:

  Panel A  "MG-state enrichment per module"
      rows = 12 MG-state / lipid-MG gene sets ; cols = WGCNA modules
      cell = overlap COUNT (bold) + BH-FDR stars; fill = log2 fold-enrichment
      (blue-white-red, ±3); significant cells (padj<0.05) black-bordered.
      Source: MG_focused_enrichment_10K_p12.csv (one-sided Fisher / ORA).

  Panel B  "WGCNA module-trait correlation"
      rows = DAB1 KD | DAB2 KD | ε4 vs ε3 (Ctrl) ; cols = SAME modules
      cell = Pearson / point-biserial r + BH-FDR stars; fill ±1.
      Source: module_trait_correlations_expanded_10K_p12.csv.

MAIN : the 7 module-trait-significant modules as columns.
SUPP : all 14 modules (grey excluded) as columns.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

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

ROOT = Path("<EDIT_PROJECT_ROOT>/Final Figures /WGCNA_10K")
ENR_CSV = ROOT / "tables" / "MG_focused_enrichment_10K_p12.csv"
MT_CSV = ROOT / "tables" / "module_trait_correlations_expanded_10K_p12.csv"
OUT_DIR = ROOT / "figures"; OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = ROOT / "figure_data"; DATA_DIR.mkdir(parents=True, exist_ok=True)
N_UNIV = 10000

MOD_COL = {
    "turquoise": "#1FB6B6", "blue": "#2C6FB3", "brown": "#7F4F2D",
    "yellow": "#C9A800", "green": "#2E9B3C", "red": "#D62728",
    "black": "#333333", "pink": "#E377A8", "magenta": "#C026B3",
    "purple": "#7E3FAE", "greenyellow": "#7FA300", "tan": "#B08A4F",
    "salmon": "#E8735E", "cyan": "#159AA8",
}
MOD_ANNO = {
    "brown": "Metabolic /\nviral carcinog.", "green": "Circadian /\ngene expr.",
    "magenta": "Innate immune\n/ complement", "blue": "Ribosome /\nmetabolic",
    "greenyellow": "Lysosome /\nlipid", "salmon": "Innate immune\n/ lysosomal",
    "purple": "Lysosome /\nDAB2·TREM2", "turquoise": "Tube morph. /\ncytoskel.",
    "yellow": "Nucleosome /\nchromatin", "pink": "Inflammatory\n/ IFN",
    "red": "Cation\ntransport", "tan": "IFN / innate\nimmune",
    "cyan": "MERTK /\nLRP5", "black": "Antiviral /\nIFN (TLR)",
}

# significant module column order (DAB2-r desc, grey excluded)
SIG_MODS = ["brown", "green", "magenta", "blue", "greenyellow", "salmon",
            "purple"]
ALL_MODS = ["brown", "green", "pink", "yellow", "black", "tan", "turquoise",
            "red", "cyan", "magenta", "blue", "greenyellow", "salmon",
            "purple"]

# Panel A rows: (gene_set key, display label, group)
ROWS_A = [
    ("MG_HM",  "Homeostatic (HM)",          "core"),
    ("MG_DAM", "Disease-associated (DAM)",   "core"),
    ("MG_CRM", "Cytokine-response (CRM)",    "core"),
    ("MG_IRM", "Interferon-response (IRM)",  "core"),
    ("MG_HLA", "Antigen-presenting (HLA)",   "core"),
    ("MG_RM",  "Ribosomal (RM)",             "core"),
    ("Trem2_LAM_core",             "TREM2-LAM core",          "lipid"),
    ("Nugent2020_TREM2_lipid_MG",  "TREM2-lipid MG (Nugent)", "lipid"),
    ("LXR_cholesterol_efflux_MG",  "LXR cholesterol-efflux",  "lipid"),
    ("Haney2024_LDAM_markers",     "LDAM (Haney)",            "lipid"),
    ("Marschallinger2020_LDAM_up", "LDAM (Marschallinger)",   "lipid"),
    ("Victor2022_APOE4_lipid_iMG", "APOE4-lipid iMG (Victor)", "lipid"),
]
GROUP_LABEL = {"core": "Mancuso 2024\ncore MG states",
               "lipid": "LAM / lipid\nMG states"}
TRAITS_B = [("DAB1 KD", "DAB1siRNA"), ("DAB2 KD", "DAB2siRNA"),
            ("ε4 vs ε3 (Ctrl)", "E4_Ctrl")]


def stars_padj(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


enr = pd.read_csv(ENR_CSV)
mt = pd.read_csv(MT_CSV)
R_tab = mt.pivot(index="module", columns="trait", values="r")
P_tab = mt.pivot(index="module", columns="trait", values="padj")

NORM_A = TwoSlopeNorm(vcenter=0.0, vmin=-3, vmax=3)
NORM_B = TwoSlopeNorm(vcenter=0.0, vmin=-1.0, vmax=1.0)


def render(modules, mode):
    n_cols = len(modules)
    nA, nB = len(ROWS_A), len(TRAITS_B)

    # ---- assemble panel A ----
    A_count = np.zeros((nA, n_cols), int)
    A_log2 = np.full((nA, n_cols), np.nan)
    A_padj = np.full((nA, n_cols), np.nan)
    A_size = np.zeros(nA, int)
    for i, (key, _, _) in enumerate(ROWS_A):
        sub = enr[enr.gene_set == key]
        if not sub.empty:
            A_size[i] = int(sub.set_size_bg.iloc[0])
        for j, mod in enumerate(modules):
            row = sub[sub.module == mod]
            if row.empty:
                continue
            row = row.iloc[0]
            A_count[i, j] = int(row.overlap)
            fe = float(row.fold_enrichment)
            A_log2[i, j] = np.log2(fe) if fe > 0 else np.nan
            A_padj[i, j] = float(row.padj)

    # ---- assemble panel B ----
    B_r = np.full((nB, n_cols), np.nan)
    B_p = np.full((nB, n_cols), np.nan)
    for i, (_, key) in enumerate(TRAITS_B):
        for j, mod in enumerate(modules):
            if mod in R_tab.index and key in R_tab.columns:
                B_r[i, j] = R_tab.loc[mod, key]
                B_p[i, j] = P_tab.loc[mod, key]

    # ---- geometry ----
    CELL_W = 1.95 if mode == "main" else 1.55
    CELL_H = 1.75
    XC = [j * CELL_W for j in range(n_cols)]
    XLIM = (-CELL_W/2 - 0.2, (n_cols-1)*CELL_W + CELL_W/2 + 0.2)

    fig_w_cm = 7.5 + n_cols * CELL_W + 3.5
    fig_h_cm = 4.5 + (nA + nB) * CELL_H + 6.0
    fig = plt.figure(figsize=(fig_w_cm * CM, fig_h_cm * CM))

    left = 6.8 / fig_w_cm
    right = 1 - 3.2 / fig_w_cm
    top = 1 - 1.9 / fig_h_cm
    bottom_B = 5.2 / fig_h_cm
    plot_w = right - left
    gap = 0.055
    usable_h = top - bottom_B - gap
    hA = usable_h * (nA / (nA + nB))
    hB = usable_h * (nB / (nA + nB))
    axA = fig.add_axes([left, top - hA, plot_w, hA])
    axB = fig.add_axes([left, bottom_B, plot_w, hB])

    def cell(ax, xc, yc, fc, ec="white", lw=1.1):
        ax.add_patch(Rectangle((xc-CELL_W/2, yc-CELL_H/2), CELL_W, CELL_H,
                     facecolor=fc, edgecolor=ec, linewidth=lw))

    # ---- Panel A ----
    cl = np.clip(A_log2, -3, 3)
    for i in range(nA):
        yc = i * CELL_H
        for j in range(n_cols):
            xc = XC[j]; v = cl[i, j]; a = A_count[i, j]; p = A_padj[i, j]
            if np.isnan(v) or a == 0:
                cell(axA, xc, yc, "white", ec="#CFCFCF", lw=0.7); continue
            cell(axA, xc, yc, NES_CMAP(NORM_A(v)))
            s = stars_padj(p)
            tcol = "white" if abs(v) > 1.8 else "#111"
            if s:
                axA.text(xc, yc - 0.13, f"{a}", ha="center", va="center",
                         fontsize=F["value"]+3, color=tcol, fontweight="bold")
                axA.text(xc, yc + 0.46, s, ha="center", va="center",
                         fontsize=F["star"], color=tcol, fontweight="bold")
                axA.add_patch(Rectangle((xc-CELL_W/2, yc-CELL_H/2), CELL_W,
                              CELL_H, facecolor="none", edgecolor="black",
                              linewidth=1.5, zorder=6))
            else:
                axA.text(xc, yc, f"{a}", ha="center", va="center",
                         fontsize=F["value"]+3, color=tcol, fontweight="bold")
    axA.set_yticks([i*CELL_H for i in range(nA)])
    axA.set_yticklabels([f"{lbl}  ({A_size[i]})" for i, (_, lbl, _) in
                         enumerate(ROWS_A)], fontsize=F["row"]-3.0)
    axA.tick_params(axis="y", left=False, pad=5)
    axA.set_xticks(XC)
    axA.set_xticklabels([m.capitalize() for m in modules],
                        fontsize=F["col"]-2.5, fontweight="bold")
    for tick, m in zip(axA.get_xticklabels(), modules):
        tick.set_color(MOD_COL[m])
    axA.tick_params(axis="x", top=True, labeltop=True, bottom=False,
                    labelbottom=False, length=0, pad=5)
    for sp in axA.spines.values():
        sp.set_visible(False)
    axA.set_xlim(*XLIM)
    axA.set_ylim((nA-1)*CELL_H + CELL_H/2 + 0.2, -CELL_H/2 - 0.2)

    # group brackets on the far left
    gtr = axA.get_yaxis_transform()
    spans = {}
    for i, (_, _, grp) in enumerate(ROWS_A):
        spans.setdefault(grp, [i, i])[1] = i
    x_br = -0.30 if mode == "main" else -0.235
    for grp, (i0, i1) in spans.items():
        y0 = i0*CELL_H - CELL_H/2; y1 = i1*CELL_H + CELL_H/2
        axA.plot([x_br, x_br], [y0, y1], transform=gtr, color="#555",
                 lw=1.4, clip_on=False)
        axA.text(x_br - 0.012, (y0+y1)/2, GROUP_LABEL[grp], transform=gtr,
                 ha="right", va="center", rotation=90, fontsize=F["col"]-2.5,
                 fontweight="bold", color="#555", linespacing=0.95)

    # ---- Panel B ----
    for i in range(nB):
        yc = i * CELL_H
        for j in range(n_cols):
            xc = XC[j]; r = B_r[i, j]; p = B_p[i, j]
            if np.isnan(r):
                cell(axB, xc, yc, "white", ec="#CFCFCF", lw=0.7); continue
            cell(axB, xc, yc, NES_CMAP(NORM_B(r)))
            s = stars_padj(p)
            tcol = "white" if abs(r) > 0.6 else "#111"
            if s:
                axB.text(xc, yc - 0.40, s, ha="center", va="center",
                         fontsize=F["star"], color=tcol, fontweight="bold")
                axB.text(xc, yc + 0.20, f"{r:+.2f}", ha="center", va="center",
                         fontsize=F["value"]+1, color=tcol, fontweight="bold")
                axB.add_patch(Rectangle((xc-CELL_W/2, yc-CELL_H/2), CELL_W,
                              CELL_H, facecolor="none", edgecolor="black",
                              linewidth=1.5, zorder=6))
            else:
                axB.text(xc, yc, f"{r:+.2f}", ha="center", va="center",
                         fontsize=F["value"]+1, color=tcol, fontweight="bold")
    axB.set_yticks([i*CELL_H for i in range(nB)])
    axB.set_yticklabels([t[0] for t in TRAITS_B], fontsize=F["row"]-2.0)
    axB.tick_params(axis="y", left=False, pad=5)
    axB.set_xticks(XC); axB.set_xticklabels([""]*n_cols)
    axB.tick_params(axis="x", length=0)
    for sp in axB.spines.values():
        sp.set_visible(False)
    axB.set_xlim(*XLIM)
    axB.set_ylim((nB-1)*CELL_H + CELL_H/2 + 0.2, -CELL_H/2 - 0.2)

    # bottom 2-line module labels
    fig.canvas.draw()
    posB = axB.get_position()
    xspan = XLIM[1] - XLIM[0]
    y_name = posB.y0 - 0.030
    y_anno = posB.y0 - 0.052
    for j, mod in enumerate(modules):
        fx = posB.x0 + (XC[j] - XLIM[0]) / xspan * posB.width
        fig.text(fx, y_name, mod, ha="center", va="center",
                 fontsize=F["row"]-3.0, fontweight="bold", color=MOD_COL[mod])
        fig.text(fx, y_anno, MOD_ANNO.get(mod, ""), ha="center", va="top",
                 fontsize=F["value"]-2.0, fontweight="bold",
                 color=MOD_COL[mod], linespacing=1.05)

    # headers
    def header(ax, text):
        pos = ax.get_position(); bh = 0.030; by = pos.y1 + 0.010
        fig.patches.append(FancyBboxPatch((pos.x0, by), pos.width, bh,
                           boxstyle="round,pad=0.004,rounding_size=0.012",
                           transform=fig.transFigure, facecolor="#EBEBEB",
                           edgecolor="none", clip_on=False, zorder=10))
        fig.text(pos.x0 + pos.width/2, by + bh/2, text, ha="center",
                 va="center", fontsize=F["title"]-2, fontweight="bold",
                 color="#222", zorder=11)
    header(axA, "A.  MG-state enrichment per module")
    header(axB, "B.  WGCNA module–trait correlation")

    # colourbars
    cbA = fig.add_axes([right + 0.018, top - hA + hA*0.12, 0.011, hA*0.55])
    c1 = mpl.colorbar.ColorbarBase(cbA, cmap=NES_CMAP, norm=NORM_A,
                                   ticks=[-3, -1.5, 0, 1.5, 3])
    c1.set_label("log2 FE", fontsize=F["cbar_lbl"]-2, fontweight="bold")
    c1.ax.tick_params(labelsize=F["cbar_tick"]-1); c1.outline.set_linewidth(0.4)
    cbB = fig.add_axes([right + 0.018, bottom_B + hB*0.10, 0.011, hB*0.7])
    c2 = mpl.colorbar.ColorbarBase(cbB, cmap=NES_CMAP, norm=NORM_B,
                                   ticks=[-1, 0, 1])
    c2.set_label("Pearson r", fontsize=F["cbar_lbl"]-2, fontweight="bold")
    c2.ax.tick_params(labelsize=F["cbar_tick"]-1); c2.outline.set_linewidth(0.4)

    ttl = "WGCNA 10K — module × MG-state map"
    fig.text(left, 0.992, ttl + ("" if mode == "main" else " (supp.)"),
             ha="left", va="top", fontsize=F["title"]-2, fontweight="bold")
    foot = (
        "Panel A: one-sided Fisher exact / ORA; cell number = overlapping "
        "genes, fill = log2 fold-enrichment (±3), black border + stars = "
        "BH-FDR sig; row label = gene-set size in 10,000-gene universe.  "
        "Panel B: point-biserial / Pearson r (DAB1 KD pooled n=25, DAB2 KD "
        "pooled n=22, ε4 vs ε3 Ctrl-only n=8).  ***padj<0.001, **<0.01, "
        "*<0.05.  "
        + ("MAIN: 7 module-trait-significant modules."
           if mode == "main" else "SUPP: all 14 modules (grey excluded)."))
    fig.text(0.012, 0.010, foot, ha="left", va="bottom",
             fontsize=F["caption"]-2.0, color="#555", style="italic", wrap=True)

    out = OUT_DIR / (f"Fig3_10K_MGstate_composite_"
                     f"{'main' if mode == 'main' else 'supp'}.png")
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"saved {out.name}  (cols={n_cols})")


if __name__ == "__main__":
    render(SIG_MODS, "main")
    render(ALL_MODS, "supp")
