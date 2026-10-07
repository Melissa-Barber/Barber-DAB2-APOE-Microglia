#!/usr/bin/env python3
"""
WGCNA top-10,000-gene run (n=39 BIONMG iPSC-MG) — microglia / AD-relevant
PATHWAY enrichment per module, single-panel horizontal fold-enrichment bars
grouped & coloured by module (matches the 5k pathway_barplot style).

Source: MG_focused_enrichment_10K_p12.csv (one-sided Fisher / ORA:
fold_enrichment, overlap, padj).  Only canonical functional pathways are used
here — the Mancuso MG_* cell-state sets and the author-named MG/lipid
signatures (Haney/Marschallinger/Nugent/TREM2-LAM/Victor/LXR) are excluded
(those belong to the module x MG-state figure).

MAIN  : the 7 module-trait-significant modules only.
SUPP  : all modules with >=1 significant functional pathway.
Per module: top 6 significant pathways (padj<0.05, overlap>=4), redundant
sets collapsed by gene Jaccard>0.6, ranked by significance.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CM = 1 / 2.54
DPI = 600
FONTS = {"title": 14, "row": 13, "col": 12, "value": 9, "caption": 9.5}
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10, "axes.linewidth": 0.6,
    "savefig.dpi": DPI, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42})
F = FONTS

ROOT = Path("<EDIT_PROJECT_ROOT>/Final Figures /WGCNA_10K")
SRC = ROOT / "tables" / "MG_focused_enrichment_10K_p12.csv"
OUT_DIR = ROOT / "figures"; OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = ROOT / "figure_data"; DATA_DIR.mkdir(parents=True, exist_ok=True)

MOD_COL = {
    "turquoise": "#1FB6B6", "blue": "#2C6FB3", "brown": "#7F4F2D",
    "yellow": "#C9A800", "green": "#2E9B3C", "red": "#D62728",
    "black": "#333333", "pink": "#E377A8", "magenta": "#C026B3",
    "purple": "#7E3FAE", "greenyellow": "#7FA300", "tan": "#B08A4F",
    "salmon": "#E8735E", "cyan": "#159AA8",
}
# legend annotation reflecting the pathways shown in THIS figure
MOD_LEG = {
    "blue":        "Metabolic / lipid (mTORC1, sterol)",
    "magenta":     "Innate immune / complement",
    "salmon":      "Innate immune / lysosomal",
    "purple":      "Lysosomal / phagocytic (DAB2)",
    "greenyellow": "Lysosomal / lipid",
    "brown":       "TNFα–NFκB / metabolic",
    "green":       "TNFα–NFκB / circadian",
    "black":       "Antiviral / IFN",
    "pink":        "Inflammatory / IFN",
    "tan":         "IFN / innate immune",
}

SIG_MODS = ["blue", "magenta", "salmon", "purple", "greenyellow",
            "brown", "green"]              # module-trait-significant (MAIN order)

STATES = {"MG_CRM", "MG_CRM1", "MG_CRM2", "MG_DAM", "MG_HLA", "MG_HM",
          "MG_IRM", "MG_RM", "MG_tCRM"}
NAMED = {"Haney2024_LDAM_markers", "LXR_cholesterol_efflux_MG",
         "Marschallinger2020_LDAM_up", "Nugent2020_TREM2_lipid_MG",
         "Trem2_LAM_core", "Victor2022_APOE4_lipid_iMG"}

LABELS = {
    "ADAPTIVE_IMMUNE_SYSTEM": "Adaptive immune system",
    "APOPTOSIS": "Apoptosis",
    "BILE_ACID_METABOLISM": "Bile-acid metabolism",
    "CHEMOKINE_SIGNALING_PATHWAY": "Chemokine signaling",
    "CHOLESTEROL_BIOSYNTHESIS": "Cholesterol biosynthesis",
    "CHOLESTEROL_HOMEOSTASIS": "Cholesterol homeostasis",
    "COAGULATION": "Coagulation",
    "COMPLEMENT": "Complement",
    "COMPLEMENT_AND_COAGULATION_CASCADES": "Complement & coagulation cascades",
    "CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION": "Cytokine–receptor interaction",
    "CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM": "Cytokine signaling (immune)",
    "FATTY_ACID_METABOLISM": "Fatty-acid metabolism",
    "FCGR_ACTIVATION": "FcγR activation",
    "FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS": "Fcγ-R–mediated phagocytosis",
    "IL6_JAK_STAT3_SIGNALING": "IL6–JAK–STAT3 signaling",
    "INFLAMMATORY_RESPONSE": "Inflammatory response",
    "INNATE_IMMUNE_SYSTEM": "Innate immune system",
    "INTERFERON_ALPHA_RESPONSE": "Interferon-α response",
    "INTERFERON_GAMMA_RESPONSE": "Interferon-γ response",
    "LYSOSOME": "Lysosome",
    "LYSOSOME_VESICLE_BIOGENESIS": "Lysosome / vesicle biogenesis",
    "MTORC1_SIGNALING": "mTORC1 signaling",
    "OXIDATIVE_PHOSPHORYLATION": "Oxidative phosphorylation",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY": "Reactive oxygen species",
    "SIGNALING_BY_INTERLEUKINS": "Signaling by interleukins",
    "TNFA_SIGNALING_VIA_NFKB": "TNFα signaling via NF-κB",
    "TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY": "Toll-like receptor signaling",
    "XENOBIOTIC_METABOLISM": "Xenobiotic metabolism",
}

PADJ_CUT, MIN_OVERLAP, JACCARD_CUT, N_TOP = 0.05, 4, 0.6, 6


def clean(gs):
    return LABELS.get(gs, gs.replace("_", " ").capitalize())


def fmt_q(q):
    return (f"{q:.0e}".replace("e-0", "e-").replace("e+0", "e+")
            if q < 0.01 else f"{q:.2f}")


def genes(s):
    return set(str(s).replace(",", ";").split(";")) if pd.notna(s) else set()


df = pd.read_csv(SRC)
func = df[~df.gene_set.isin(STATES | NAMED)].copy()


def select(modules):
    chosen = []
    for mod in modules:
        sub = (func[(func.module == mod) & (func.padj < PADJ_CUT)
                    & (func.overlap >= MIN_OVERLAP)]
               .sort_values("padj"))
        kept, kept_g = [], []
        for _, r in sub.iterrows():
            g = genes(r["overlap_genes"])
            if any(len(g & ks) / max(1, len(g | ks)) > JACCARD_CUT
                   for ks in kept_g):
                continue
            kept.append(r); kept_g.append(g)
            if len(kept) >= N_TOP:
                break
        if kept:
            chosen.append(pd.DataFrame(kept))
    return pd.concat(chosen, ignore_index=True) if chosen else pd.DataFrame()


def render(modules, title, footer, out_name, data_name):
    sel = select(modules)
    rows, group_spans, GAP, y = [], [], 0.55, 0.0
    used_mods = []
    for mod in modules:
        sub = sel[sel.module == mod].sort_values("padj")
        if sub.empty:
            continue
        used_mods.append(mod)
        for _, r in sub.iterrows():
            rows.append(dict(y=y, module=mod, label=clean(r["gene_set"]),
                             fold=r["fold_enrichment"], count=int(r["overlap"]),
                             q=r["padj"]))
            y += 1.0
        y += GAP
    total_y = y
    fold_max = max(r["fold"] for r in rows)
    x_max = fold_max * 1.30

    fig_w_cm = 30.0
    fig_h_cm = 3.0 + 0.55 * total_y
    fig = plt.figure(figsize=(fig_w_cm * CM, fig_h_cm * CM))
    gs = fig.add_gridspec(1, 3, width_ratios=[11.0, 13.0, 6.0], wspace=0.03,
                          left=0.012, right=0.988,
                          top=1 - 1.7 / fig_h_cm, bottom=2.1 / fig_h_cm)

    ax_lbl = fig.add_subplot(gs[0])
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(-0.5 - GAP, total_y - 0.5)
    ax_lbl.invert_yaxis(); ax_lbl.axis("off")
    for r in rows:
        ax_lbl.text(0.985, r["y"], r["label"], ha="right", va="center",
                    fontsize=F["row"] - 3.5, color=MOD_COL[r["module"]],
                    weight="bold")

    ax_bar = fig.add_subplot(gs[1], sharey=ax_lbl)
    for r in rows:
        ax_bar.barh(r["y"], r["fold"], height=0.78, color=MOD_COL[r["module"]],
                    edgecolor="#333", linewidth=0.4, zorder=3)
    ax_bar.set_xlim(0, x_max); ax_bar.set_ylim(-0.5 - GAP, total_y - 0.5)
    ax_bar.invert_yaxis(); ax_bar.set_yticks([])
    for sp in ("top", "right", "left"):
        ax_bar.spines[sp].set_visible(False)
    ax_bar.tick_params(axis="x", labelsize=F["col"] - 1)
    ax_bar.set_xlabel("Fold enrichment", fontsize=F["col"], labelpad=9)
    ax_bar.grid(axis="x", color="#dddddd", linewidth=0.5, zorder=0)
    off = fold_max * 0.02
    for r in rows:
        ax_bar.text(r["fold"] + off, r["y"], f"n={r['count']}  q={fmt_q(r['q'])}",
                    ha="left", va="center", fontsize=F["value"] - 1.0,
                    color="#222", zorder=4)

    ax_leg = fig.add_subplot(gs[2])
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1); ax_leg.axis("off")
    ax_leg.text(0.0, 0.985, "WGCNA module", ha="left", va="top",
                fontsize=F["col"], weight="bold", color="#222")
    sw_h, y0, dy = 0.045, 0.90, 0.085
    for i, mod in enumerate(used_mods):
        yy = y0 - i * dy
        ax_leg.add_patch(Rectangle((0.0, yy - sw_h / 2), 0.055, sw_h,
                         facecolor=MOD_COL[mod], edgecolor="#333",
                         linewidth=0.4, transform=ax_leg.transAxes))
        ax_leg.text(0.085, yy, f"{mod} ({MOD_LEG.get(mod, '')})", ha="left",
                    va="center", fontsize=F["value"] - 0.5,
                    color=MOD_COL[mod], weight="bold")

    fig.suptitle(title, x=0.012, y=0.992, ha="left",
                 fontsize=F["title"], weight="bold")
    fig.text(0.012, 0.008, footer, ha="left", va="bottom",
             fontsize=F["caption"] - 2.0, color="#555", style="italic")
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    sel.assign().to_csv(DATA_DIR / data_name, index=False)
    print(f"saved {out.name}  ({len(rows)} bars; modules={used_mods})")


FOOT = ("One-sided Fisher exact / ORA, BH-FDR.  Bar = fold enrichment; "
        "n = module genes in pathway (overlap); q = BH-adjusted p.  Top 6 "
        "significant pathways per module (padj<0.05, overlap≥4), redundant "
        "sets collapsed (gene Jaccard>0.6).  Universe = 10,000 WGCNA-expressed "
        "genes.  Mancuso MG-state & author-named lipid-MG signatures excluded "
        "(shown in the module × MG-state figure).  n=39 BIONMG iPSC-MG.")

if __name__ == "__main__":
    render(SIG_MODS,
           "Microglia / AD-relevant pathway enrichment — significant WGCNA "
           "modules (top-10,000-gene run, n=39)",
           FOOT + "  MAIN: 7 module-trait-significant modules.",
           "Fig2_10K_pathway_barplot_main.png",
           "Fig2_10K_pathway_bars_main.csv")
    all_mods = ["blue", "magenta", "salmon", "purple", "greenyellow",
                "brown", "green", "black", "pink", "tan"]
    render(all_mods,
           "Microglia / AD-relevant pathway enrichment — all WGCNA modules "
           "(top-10,000-gene run, n=39)",
           FOOT + "  SUPP: all modules with ≥1 significant functional pathway.",
           "Fig2_10K_pathway_barplot_supp.png",
           "Fig2_10K_pathway_bars_supp.csv")
