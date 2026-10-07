#!/usr/bin/env python3
"""
WGCNA top-10,000-gene run (n=39 BIONMG iPSC-MG) — pathway enrichment of the
modules INDUCED vs SUPPRESSED by DAB2 knockdown.

Companion to the DAB2-KD module-eigengene correlation bar plot: the 7
module-trait-significant modules are split by the sign of their point-biserial
correlation with DAB2 KD (Fig_DAB2KD_module_switch_10K.csv) into

    Induced  by DAB2 KD (module ↑) : brown (+0.80***), green (+0.56*)
    Suppressed by DAB2 KD (module ↓): purple (-0.83***), salmon (-0.82***),
                                       greenyellow (-0.67**), blue (-0.49*),
                                       magenta (-0.47*)

and within each group the per-module functional pathway enrichment is shown as
horizontal fold-enrichment bars, coloured by module. Each module is annotated
with its DAB2-KD correlation r and significance.

TWO versions:
  FOCUSED  (main) — curated microglia / lipid / AD functional gene-set panel
                    (MG_focused_enrichment_10K_p12.csv; Mancuso MG-state and
                    author-named lipid-MG signatures excluded).
  UNBIASED (supp) — database-wide GO-BP + KEGG ORA, no keyword filter
                    (GO_BP_enrichment_10K_p12.csv + KEGG_enrichment_10K_p12.csv).

Redundancy collapse: the manuscript-standard greedy Jaccard > 0.4, applied — as
for every WGCNA ORA panel — to the genes driving each term within its module
(the module-overlap genes), keeping the most-significant representative.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch

CM = 1 / 2.54
DPI = 600
FONTS = {"title": 14, "group": 13, "mod": 11.5, "row": 12, "value": 9,
         "caption": 9.5, "leg": 9}
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10, "axes.linewidth": 0.6,
    "savefig.dpi": DPI, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42})
F = FONTS

ROOT = Path("<EDIT_PROJECT_ROOT>/Final Figures /WGCNA_10K")
TAB = ROOT / "tables"
OUT_DIR = ROOT / "figures"; OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = ROOT / "figure_data"; DATA_DIR.mkdir(parents=True, exist_ok=True)

MOD_COL = {
    "turquoise": "#1FB6B6", "blue": "#2C6FB3", "brown": "#7F4F2D",
    "yellow": "#C9A800", "green": "#2E9B3C", "red": "#D62728",
    "black": "#333333", "pink": "#E377A8", "magenta": "#C026B3",
    "purple": "#7E3FAE", "greenyellow": "#7FA300", "tan": "#B08A4F",
    "salmon": "#E8735E", "cyan": "#159AA8",
}
MOD_LEG = {
    "blue":        "Metabolic / lipid (mTORC1, sterol)",
    "magenta":     "Innate immune / complement",
    "salmon":      "Innate immune / lysosomal (HM)",
    "purple":      "Lysosomal / phagocytic (DAM·TREM2)",
    "greenyellow": "Lysosomal / lipid (DAM-assoc.)",
    "brown":       "TNFα–NFκB / metabolic (CRM)",
    "green":       "TNFα–NFκB / circadian",
}

# DAB2-KD module-eigengene point-biserial r + BH-FDR stars, read straight from
# the source table so the displayed values match exactly (no transcription).
_dab2 = pd.read_csv(ROOT / "figure_data" / "Fig_DAB2KD_module_switch_10K.csv")
DAB2 = {row.module: (float(row.r), str(row.sig))
        for row in _dab2.itertuples()}
# induced = r>0 (desc); suppressed = r<0 (most-negative first)
INDUCED    = [m for m in _dab2.sort_values("r", ascending=False).module
              if DAB2[m][0] > 0]
SUPPRESSED = [m for m in _dab2.sort_values("r", ascending=True).module
              if DAB2[m][0] < 0]
GROUPS = [("Induced by DAB2 knockdown   (module ↑)", INDUCED),
          ("Suppressed by DAB2 knockdown   (module ↓)", SUPPRESSED)]

# Mancuso MG-state sets & author-named lipid-MG signatures — excluded from the
# functional focused panel (they belong to the module × MG-state figure).
STATES = {"MG_CRM", "MG_CRM1", "MG_CRM2", "MG_DAM", "MG_HLA", "MG_HM",
          "MG_IRM", "MG_RM", "MG_tCRM"}
NAMED = {"Haney2024_LDAM_markers", "LXR_cholesterol_efflux_MG",
         "Marschallinger2020_LDAM_up", "Nugent2020_TREM2_lipid_MG",
         "Trem2_LAM_core", "Victor2022_APOE4_lipid_iMG"}

LABELS = {
    "ADAPTIVE_IMMUNE_SYSTEM": "Adaptive immune system",
    "APOPTOSIS": "Apoptosis", "BILE_ACID_METABOLISM": "Bile-acid metabolism",
    "CHEMOKINE_SIGNALING_PATHWAY": "Chemokine signaling",
    "CHOLESTEROL_BIOSYNTHESIS": "Cholesterol biosynthesis",
    "CHOLESTEROL_HOMEOSTASIS": "Cholesterol homeostasis",
    "COAGULATION": "Coagulation", "COMPLEMENT": "Complement",
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

PADJ_CUT, MIN_OVERLAP, JACCARD_CUT, N_TOP = 0.05, 4, 0.4, 5


def fmt_q(q):
    return (f"{q:.0e}".replace("e-0", "e-").replace("e+0", "e+")
            if q < 0.01 else f"{q:.2f}")


def genes(s):
    return set(str(s).replace(",", ";").replace("/", ";").split(";")) \
        if pd.notna(s) else set()


def clean_focused(gs):
    return LABELS.get(gs, gs.replace("_", " ").capitalize())


def clean_unbiased(desc):
    s = str(desc).replace("_", " ").strip()
    if s.isupper() or s.islower():
        s = s[:1].upper() + s[1:]
    s = re.sub(r"\bnf kappa b\b", "NF-κB", s, flags=re.I)
    s = re.sub(r"\btnf\b", "TNF", s, flags=re.I)
    s = re.sub(r"\bmhc\b", "MHC", s, flags=re.I)
    return s


# ------------------------------------------------------------------ loaders
def load_focused():
    df = pd.read_csv(TAB / "MG_focused_enrichment_10K_p12.csv")
    df = df[~df.gene_set.isin(STATES | NAMED)].copy()
    df = df.rename(columns={"gene_set": "term", "fold_enrichment": "fold",
                            "overlap": "count", "overlap_genes": "ov"})
    df["clean"] = df["term"].apply(clean_focused)
    return df[["module", "term", "clean", "fold", "count", "padj", "ov"]]


def load_unbiased():
    frames = []
    go = pd.read_csv(TAB / "GO_BP_enrichment_10K_p12.csv")
    go = go.rename(columns={"Description": "term", "FoldEnrichment": "fold",
                            "Count": "count", "p.adjust": "padj",
                            "geneID": "ov"})
    go["src"] = "GO:BP"
    frames.append(go[["module", "term", "fold", "count", "padj", "ov", "src"]])
    kg = pd.read_csv(TAB / "KEGG_enrichment_10K_p12.csv")
    kg = kg.rename(columns={"Description": "term", "FoldEnrichment": "fold",
                            "Count": "count", "p.adjust": "padj",
                            "geneID": "ov"})
    kg["src"] = "KEGG"
    frames.append(kg[["module", "term", "fold", "count", "padj", "ov", "src"]])
    df = pd.concat(frames, ignore_index=True)
    df["clean"] = df["term"].apply(clean_unbiased)
    return df


# ------------------------------------------------------------- term selection
def select(df, modules):
    """Per module: sig terms (padj<cut, count>=MIN_OVERLAP), greedy overlap-gene
    Jaccard>0.4 collapse (most-significant rep kept), top N_TOP."""
    out = {}
    for mod in modules:
        sub = (df[(df.module == mod) & (df.padj < PADJ_CUT)
                  & (df["count"] >= MIN_OVERLAP)]
               .sort_values("padj")
               .drop_duplicates(subset="clean", keep="first"))
        kept, kept_g = [], []
        for _, r in sub.iterrows():
            g = genes(r["ov"])
            if any(len(g & ks) / max(1, len(g | ks)) > JACCARD_CUT
                   for ks in kept_g):
                continue
            kept.append(r); kept_g.append(g)
            if len(kept) >= N_TOP:
                break
        if kept:
            out[mod] = pd.DataFrame(kept)
    return out


# -------------------------------------------------------------------- render
def render(sel_by_mod, title, footer, out_name, data_name):
    # Build the item list (group headers, module headers, bars) top -> bottom
    items, y = [], 0.0
    MOD_GAP, GRP_GAP, GRP_H, MOD_H, BAR_H = 0.45, 1.05, 1.15, 1.05, 1.0
    export = []
    for gname, gmods in GROUPS:
        present = [m for m in gmods if m in sel_by_mod and len(sel_by_mod[m])]
        if not present:
            continue
        items.append(dict(t="grp", y=y, text=gname)); y += GRP_H
        for m in present:
            r, sig = DAB2[m]
            items.append(dict(t="mod", y=y, module=m,
                              text=f"{m}  —  {MOD_LEG.get(m, '')}",
                              rtext=f"DAB2 r = {r:+.2f} {sig}")); y += MOD_H
            for _, row in sel_by_mod[m].sort_values("padj").iterrows():
                items.append(dict(t="bar", y=y, module=m, label=row["clean"],
                                  fold=row["fold"], count=int(row["count"]),
                                  q=row["padj"])); y += BAR_H
                export.append(dict(group=gname, module=m, DAB2_r=r, DAB2_sig=sig,
                                   pathway=row["clean"], fold_enrichment=row["fold"],
                                   overlap=int(row["count"]), padj=row["padj"]))
            y += MOD_GAP
        y += GRP_GAP
    total_y = y
    bars = [it for it in items if it["t"] == "bar"]
    fold_max = max(b["fold"] for b in bars)
    x_max = fold_max * 1.32

    fig_w_cm = 30.0
    fig_h_cm = 3.4 + 0.52 * total_y
    fig = plt.figure(figsize=(fig_w_cm * CM, fig_h_cm * CM))
    gs = fig.add_gridspec(1, 3, width_ratios=[12.5, 12.0, 5.5], wspace=0.03,
                          left=0.012, right=0.988,
                          top=1 - 1.9 / fig_h_cm, bottom=2.2 / fig_h_cm)

    ax_lbl = fig.add_subplot(gs[0])
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(-0.6, total_y - 0.4)
    ax_lbl.invert_yaxis(); ax_lbl.axis("off")

    ax_bar = fig.add_subplot(gs[1], sharey=ax_lbl)
    ax_bar.set_xlim(0, x_max); ax_bar.set_ylim(-0.6, total_y - 0.4)
    ax_bar.invert_yaxis()

    # group bands (light grey) spanning both label + bar axes
    for it in items:
        if it["t"] != "grp":
            continue
        # find the group's vertical extent
        gi = items.index(it)
        gstart = it["y"] - 0.55
        gend = total_y
        for j in range(gi + 1, len(items)):
            if items[j]["t"] == "grp":
                gend = items[j]["y"] - 0.55 - GRP_GAP
                break
        for ax in (ax_lbl, ax_bar):
            ax.axhspan(gstart, it["y"] + 0.45, color="#EFEFEF", zorder=0)

    # draw items
    for it in items:
        if it["t"] == "grp":
            ax_lbl.text(0.0, it["y"], it["text"], ha="left", va="center",
                        fontsize=F["group"], weight="bold", color="#111")
        elif it["t"] == "mod":
            col = MOD_COL[it["module"]]
            ax_lbl.add_patch(Rectangle((0.012, it["y"] - 0.30), 0.022, 0.60,
                             facecolor=col, edgecolor="#333", linewidth=0.4))
            ax_lbl.text(0.05, it["y"], it["text"], ha="left", va="center",
                        fontsize=F["mod"], weight="bold", color=col)
            ax_bar.text(x_max * 0.995, it["y"], it["rtext"], ha="right",
                        va="center", fontsize=F["mod"] - 0.5, weight="bold",
                        color=col)
        else:  # bar
            col = MOD_COL[it["module"]]
            ax_lbl.text(0.975, it["y"], it["label"], ha="right", va="center",
                        fontsize=F["row"] - 3.0, color=col)
            ax_bar.barh(it["y"], it["fold"], height=0.74, color=col,
                        edgecolor="#333", linewidth=0.4, zorder=3)
            ax_bar.text(it["fold"] + fold_max * 0.02, it["y"],
                        f"n={it['count']}  q={fmt_q(it['q'])}", ha="left",
                        va="center", fontsize=F["value"] - 1.0, color="#222",
                        zorder=4)

    ax_bar.set_yticks([])
    for sp in ("top", "right", "left"):
        ax_bar.spines[sp].set_visible(False)
    ax_bar.tick_params(axis="x", labelsize=F["value"])
    ax_bar.set_xlabel("Fold enrichment", fontsize=11, labelpad=8)
    ax_bar.grid(axis="x", color="#dddddd", linewidth=0.5, zorder=0)

    # legend column: module swatch + DAB2 r
    ax_leg = fig.add_subplot(gs[2]); ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
    ax_leg.axis("off")
    ax_leg.text(0.0, 0.99, "WGCNA module · DAB2-KD r", ha="left", va="top",
                fontsize=11, weight="bold", color="#222")
    used = [m for _, gm in GROUPS for m in gm
            if m in sel_by_mod and len(sel_by_mod[m])]
    y0, dy, sw = 0.92, 0.072, 0.05
    for i, m in enumerate(used):
        yy = y0 - i * dy
        r, sig = DAB2[m]
        ax_leg.add_patch(Rectangle((0.0, yy - sw / 2), 0.05, sw,
                         facecolor=MOD_COL[m], edgecolor="#333", linewidth=0.4))
        arrow = "↑" if r > 0 else "↓"
        ax_leg.text(0.075, yy, f"{m}  ({r:+.2f}{sig} {arrow})", ha="left",
                    va="center", fontsize=F["leg"], weight="bold",
                    color=MOD_COL[m])

    fig.suptitle(title, x=0.012, y=0.992, ha="left", fontsize=F["title"],
                 weight="bold")
    fig.text(0.012, 0.008, footer, ha="left", va="bottom",
             fontsize=F["caption"] - 2.0, color="#555", style="italic")
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    pd.DataFrame(export).to_csv(DATA_DIR / data_name, index=False)
    print(f"saved {out.name}  ({len(bars)} bars; modules={used})")


FOOT_COMMON = (
    "Modules split by the sign of their module-eigengene point-biserial "
    "correlation with DAB2 knockdown (Fig_DAB2KD_module_switch; ***padj<0.001, "
    "**<0.01, *<0.05).  Bar = fold enrichment; n = module genes in pathway "
    "(overlap); q = BH-adjusted p.  One-sided Fisher/ORA.  Top 5 significant "
    "pathways per module (padj<0.05, overlap≥4); redundant terms collapsed by "
    "overlap-gene Jaccard>0.4 (most-significant representative kept).  "
    "Universe = 10,000 WGCNA-expressed genes.  n=39 BIONMG iPSC-MG.")

if __name__ == "__main__":
    foc = select(load_focused(),
                 INDUCED + SUPPRESSED)
    render(foc,
           "Pathway enrichment of WGCNA modules induced vs suppressed by DAB2 "
           "knockdown — microglia / lipid / AD-focused",
           FOOT_COMMON + "  FOCUSED: curated MG / lipid / AD functional panel; "
           "Mancuso MG-state & author-named lipid-MG signatures excluded.",
           "Fig_10K_DAB2dir_pathway_barplot_focused_main.png",
           "Fig_10K_DAB2dir_pathway_bars_focused_main.csv")

    unb = select(load_unbiased(), INDUCED + SUPPRESSED)
    render(unb,
           "Pathway enrichment of WGCNA modules induced vs suppressed by DAB2 "
           "knockdown — unbiased (GO-BP + KEGG)",
           FOOT_COMMON + "  UNBIASED: database-wide GO-BP + KEGG ORA, no keyword "
           "filter.",
           "Fig_10K_DAB2dir_pathway_barplot_unbiased_supp.png",
           "Fig_10K_DAB2dir_pathway_bars_unbiased_supp.csv")
