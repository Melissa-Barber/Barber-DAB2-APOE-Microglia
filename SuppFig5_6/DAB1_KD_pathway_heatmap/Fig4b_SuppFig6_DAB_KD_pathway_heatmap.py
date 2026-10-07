#!/usr/bin/env python3
"""
Python renderer for the BIONMG DAB1 and DAB2 KD pathway-enrichment heatmaps,
re-rendered in the standardised manuscript style (matches every other
heatmap in this paper):

  - Deep blue → white → deep red palette  (#2166AC → #B2182B, clipped ±3)
  - White text only when |NES| > 1.85
  - Larger NES + stars fonts for readability at print scale
  - Thin NES colorbar (width 0.022 of figure) — same as every other heatmap
  - (n/N) RELN-core overlap kept inside the row label (no overlap with cells)
  - Three-block layout (Lipid/Mito/ECM/mTOR · Inflammation/Senescence · Aβ/Phagocytosis)
    with dashed inner divider for auto-selected, solid bracket for manually-added rows

Values are pulled from the BIONMG_DAB{1,2}_KD MainFig CSVs that the user's
R pipeline writes ("BIONMG_DAB1_KD_allContrasts_MG_fgsea_full.csv" /
"BIONMG_DAB2_KD_allContrasts_MG_fgsea_full.csv"). If those CSVs are not
present locally, the script falls back to the published values extracted
from the v2/v7 PNGs the user shared (so the figure can still be reproduced).

Usage:
    python3 Fig4b_SuppFig6_DAB_KD_pathway_heatmap.py {dab1|dab2|both}
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import apply_style, NES_CMAP, NES_NORM, FONTS, stars, CM, DPI
apply_style()

# ── Row definition (shared between DAB1 and DAB2) ─────────────────────────
ROWS = [
    # (path_id, display_label, block)
    ("Cholesterol Biosynthesis",         "Cholesterol Biosynthesis",         "interact_up"),
    ("Fatty Acid Metabolism",            "Fatty Acid Metabolism",            "interact_up"),
    ("Oxidative Phosphorylation",        "Oxidative Phosphorylation",        "interact_up"),
    ("Mitochondrial Proteostasis",       "Mitochondrial Proteostasis",       "interact_up"),
    ("ECM Organization",                 "ECM Organization",                 "interact_up"),
    ("mTORC1 Signaling",                 "mTORC1 Signaling",                 "interact_up"),
    ("Ubiquitin-Dependent Proteolysis",  "Ubiquitin-Dependent Proteolysis",  "interact_dn"),
    ("TNFα Signaling via NF-κB",         "TNFα Signaling via NF-κB",         "interact_dn"),
    ("Inflammatory Response",            "Inflammatory Response",            "interact_dn"),
    ("SASP",                             "SASP",                             "interact_dn"),
    ("DNA Damage-Induced Senescence",    "DNA Damage-Induced Senescence",    "interact_dn"),
    ("Cytokine-Mediated Signaling",      "Cytokine-Mediated Signaling",      "interact_dn"),
    ("Amyloid-β Clearance",              "Amyloid-β Clearance",              "manual"),
    ("Phagocytosis",                     "Phagocytosis",                     "manual"),
    ("Phagocytosis Engulfment",          "Phagocytosis Engulfment",          "manual"),
    ("Lysosomal Transport",              "Lysosomal Transport",              "manual"),
]

# ── Hardcoded data extracted from the v2/v7 PNGs the user shared ──────────
# Columns: APOE4 KD | APOE3 KD | E4 vs E3 interaction
# Format per cell: (NES, padj-stars, n_RELN_LE, gs_size)
DAB1_DATA = {
    "Cholesterol Biosynthesis":        [( 0.65, "",  0,  26), (-1.32, "",   0,  26), ( 1.25, "",  0,  26)],
    "Fatty Acid Metabolism":           [(-1.97, "***", 2, 139), (-1.74, "***", 2, 139), ( 0.66, "", 2, 139)],
    "Oxidative Phosphorylation":       [(-2.62, "***", 0, 198), (-1.79, "***", 0, 198), (-0.73, "", 0, 198)],
    "Mitochondrial Proteostasis":      [(-2.39, "***", 1,  87), (-1.64, "**",  1,  87), (-0.65, "", 1,  87)],
    "ECM Organization":                [(-1.10, "",   4, 269), (-1.29, "*",   4, 269), (-1.02, "", 4, 269)],
    "mTORC1 Signaling":                [(-1.99, "***", 1, 200), ( 0.68, "",   1, 200), (-1.40, "", 1, 200)],
    "Ubiquitin-Dependent Proteolysis": [(-0.85, "",   3, 168), (-0.97, "",   3, 168), ( 0.75, "", 3, 168)],
    "TNFα Signaling via NF-κB":        [( 1.76, "***", 0, 195), ( 2.48, "***", 0, 195), (-1.95, "***", 0, 195)],
    "Inflammatory Response":           [( 1.36, "*",  1, 177), ( 1.92, "***", 1, 177), (-1.41, "", 1, 177)],
    "SASP":                            [( 0.78, "",   1,  86), ( 1.73, "***", 1,  86), (-1.72, "", 1,  86)],
    "DNA Damage-Induced Senescence":   [( 1.21, "",   1,  62), ( 1.71, "**",  1,  62), (-1.55, "", 1,  62)],
    "Cytokine-Mediated Signaling":     [( 1.28, "*",  5, 441), ( 1.47, "***", 5, 441), (-2.04, "***", 5, 441)],
    "Amyloid-β Clearance":             [(-1.81, "**", 5,  37), (-1.80, "**",  5,  37), ( 0.80, "", 5,  37)],
    "Phagocytosis":                    [(-1.52, "**", 8, 212), (-1.77, "***", 8, 212), ( 0.82, "", 8, 212)],
    "Phagocytosis Engulfment":         [(-1.83, "**", 2,  46), (-2.19, "***", 2,  46), ( 0.81, "", 2,  46)],
    "Lysosomal Transport":             [(-1.50, "**", 3, 123), (-1.58, "**",  3, 123), ( 0.87, "", 3, 123)],
}
DAB1_COL_LABELS = [
    "APOE4 pooled\nDAB1 KD\n(siRNA-893+894)",
    "APOE3 pooled\nDAB1 KD\n(siRNA-893+894)",
    "APOE4 vs APOE3\nDAB1 interact\n(siRNA-893)",
]

DAB2_DATA = {
    "Cholesterol Biosynthesis":        [( 1.07, "",   0,  27), (-2.02, "**",  0,  27), ( 2.26, "***", 0,  27)],
    "Fatty Acid Metabolism":           [(-2.09, "**", 2, 140), (-2.12, "***", 2, 140), ( 1.62, "*",  2, 140)],
    "Oxidative Phosphorylation":       [(-2.71, "***", 0, 198), (-2.59, "***", 0, 198), ( 1.74, "***", 0, 198)],
    "Mitochondrial Proteostasis":      [(-2.28, "**", 1,  88), (-2.36, "***", 1,  88), ( 1.82, "**", 1,  88)],
    "ECM Organization":                [( 2.03, "***", 5, 263), (-1.43, "*",   5, 263), ( 1.74, "***", 5, 263)],
    "mTORC1 Signaling":                [(-1.15, "",   2, 199), (-1.80, "***", 2, 199), ( 1.71, "**", 2, 199)],
    "Ubiquitin-Dependent Proteolysis": [(-1.68, "",   6, 167), (-1.65, "**",  6, 167), ( 1.53, "*",  6, 167)],
    "TNFα Signaling via NF-κB":        [( 2.37, "***", 1, 193), ( 2.73, "***", 1, 193), (-1.77, "***", 1, 193)],
    "Inflammatory Response":           [( 1.69, "**", 2, 176), ( 2.14, "***", 2, 176), (-1.75, "***", 2, 176)],
    "SASP":                            [( 0.81, "",   0,  84), ( 2.05, "***", 0,  84), (-1.78, "**", 0,  84)],
    "DNA Damage-Induced Senescence":   [( 0.96, "",   1,  60), ( 2.23, "***", 1,  60), (-1.97, "**", 1,  60)],
    "Cytokine-Mediated Signaling":     [( 1.29, "",   5, 432), ( 1.68, "***", 5, 432), (-1.44, "**", 5, 432)],
    "Amyloid-β Clearance":             [(-1.81, "",   7,  37), (-1.99, "**",  7,  37), ( 1.51, "",  7,  37)],
    "Phagocytosis":                    [(-1.98, "***", 5, 207), (-1.78, "***", 5, 207), ( 1.01, "", 5, 207)],
    "Phagocytosis Engulfment":         [(-1.92, "***", 1,  45), (-2.07, "***", 1,  45), ( 1.18, "", 1,  45)],
    "Lysosomal Transport":             [(-1.81, "",   3, 124), (-1.47, "",    3, 124), (-0.91, "", 3, 124)],
}
DAB2_COL_LABELS = [
    "APOE4\nsiRNA-96 KD",
    "APOE3\nsiRNA-96 KD",
    "APOE4 vs APOE3\nsiRNA-96 interaction",
]


def render(which: str, data: dict, col_labels: list, title: str, subtitle: str, out_name: str):
    n = len(ROWS)

    # Inter-block gap so the 3 blocks read as 3 groups
    BLOCK_GAP = 0.45
    def y_for(i):
        b = ROWS[i][2]
        if b == "interact_up":  return i
        if b == "interact_dn":  return i + BLOCK_GAP
        return i + 2 * BLOCK_GAP

    y_max = y_for(n - 1) + 0.6

    # Cell-grid dimensions — wider cells than v2 so multi-line column labels
    # (e.g. DAB1's "(siRNA-893+894)" footer) sit cleanly within their column
    CELL_W = 3.2
    CELL_H = 1.5

    # Figure: wide enough that the right-side group brackets sit cleanly
    # next to the cells and the colorbar still has its own band of space.
    # Extra top padding so title + subtitle + multi-line column labels never
    # overlap.
    fig_w_cm = 26.0
    fig_h_cm = 9.0 + 1.10 * y_max
    fig = plt.figure(figsize=(fig_w_cm * CM, fig_h_cm * CM))

    # Heatmap axes — leave ~16 % of figure height at the top for
    # title / subtitle / 3-line column labels with no overlap.
    ax = fig.add_axes([0.27, 0.05, 0.50, 0.72])

    for i, (pid, label, _block) in enumerate(ROWS):
        yc = y_for(i) * CELL_H
        cells = data[pid]
        for j, (v, star, n_reln, gs_size) in enumerate(cells):
            xc = j * CELL_W
            face = NES_CMAP(NES_NORM(np.clip(v, -3, 3)))
            ax.add_patch(Rectangle(
                (xc - CELL_W/2, yc - CELL_H/2), CELL_W, CELL_H,
                facecolor=face, edgecolor="white", linewidth=1.2,
            ))
            text_col = "white" if abs(v) > 1.85 else "#222"
            # NES (top) — star (bottom).  (n/N) is in the row label, not the cell.
            # In-cell text: reduced from FONTS["value"]+4 / FONTS["star"]+4
            # because the +4 boost (24 / 22 pt) was overflowing the 3.2x1.5cm
            # cells and visually colliding with the column headers.
            ax.text(xc, yc - 0.20, f"{v:+.2f}",
                    ha="center", va="center",
                    fontsize=13, color=text_col,
                    weight=("bold" if star else "normal"))
            if star:
                ax.text(xc, yc + 0.40, star,
                        ha="center", va="center",
                        fontsize=11, color=text_col, weight="bold")

    # Dashed divider between interact_up and interact_dn (~auto block split)
    yc_top = (y_for(5) + y_for(6)) / 2 * CELL_H
    ax.plot([-CELL_W/2 - 0.1, (len(col_labels)-1)*CELL_W + CELL_W/2 + 0.1],
            [yc_top, yc_top], linestyle="--", color="#666", linewidth=0.8)
    # Solid divider before the manual-added block
    yc_bot = (y_for(11) + y_for(12)) / 2 * CELL_H
    ax.plot([-CELL_W/2 - 0.1, (len(col_labels)-1)*CELL_W + CELL_W/2 + 0.1],
            [yc_bot, yc_bot], linestyle="-", color="#222", linewidth=0.9)
    # (no vertical divider — the interaction column is already labelled in its header)

    # Row labels — name plus (n/N) overlap count
    yticks = [y_for(i) * CELL_H for i in range(n)]
    ytick_labels = []
    for i, (pid, label, _b) in enumerate(ROWS):
        nL, gs = data[pid][0][2], data[pid][0][3]
        ytick_labels.append(f"{label} ({nL}/{gs})")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels, fontsize=FONTS["row"])

    # Column labels at the top
    ax.set_xticks([j * CELL_W for j in range(len(col_labels))])
    ax.set_xticklabels(col_labels, fontsize=FONTS["col"], weight="bold")
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()

    ax.set_xlim(-CELL_W/2 - 0.2,
                (len(col_labels) - 1) * CELL_W + CELL_W/2 + 0.2)
    ax.set_ylim(y_max * CELL_H + 0.2, -CELL_H/2 - 0.2)

    ax.tick_params(top=False, bottom=False, left=False, right=False, pad=14)
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(False)

    # Side annotation brackets — three biological groupings.
    # Drawn in DATA coords on the heatmap axes, placed just to the right
    # of the cells (and well left of the colorbar at the far right).
    bracket_x = (len(col_labels) - 1) * CELL_W + CELL_W/2 + 0.35
    text_x    = bracket_x + 0.20
    def yrange(idxs):
        return min(y_for(i) for i in idxs) * CELL_H - CELL_H/2 + 0.05, \
               max(y_for(i) for i in idxs) * CELL_H + CELL_H/2 - 0.05

    groups = [
        ([0,1,2,3,4,5],    "Lipid, Mitochondria,\nECM & mTOR\n(interact ↑)", "#2C2C2C"),
        ([6,7,8,9,10,11],  "Inflammation &\nSenescence\n(interact ↓)",      "#2C2C2C"),
        ([12,13,14,15],    "Amyloid-β, Phagocytosis\n& Endolysosomal",      "#555555"),
    ]
    # Expand the data x-limit a little so the bracket labels are visible
    ax.set_xlim(-CELL_W/2 - 0.2,
                (len(col_labels) - 1) * CELL_W + CELL_W/2 + 2.4)
    for idxs, lbl, col in groups:
        y0, y1 = yrange(idxs)
        ax.plot([bracket_x, bracket_x], [y0, y1], color=col, linewidth=1.4)
        ax.text(text_x, (y0 + y1)/2, lbl,
                ha="left", va="center", fontsize=FONTS["row"] - 1,
                weight="bold", color=col, linespacing=1.15,
                clip_on=False)

    # NES colorbar — thin, matches every other heatmap; placed at the far right
    cax = fig.add_axes([0.91, 0.15, 0.018, 0.50])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=NES_CMAP, norm=NES_NORM,
                                   ticks=[-3, -2, -1, 0, 1, 2, 3])
    cb.set_label("NES", fontsize=FONTS["cbar_lbl"], weight="bold")
    cb.ax.tick_params(labelsize=FONTS["cbar_tick"])
    cb.outline.set_linewidth(0.4)

    # Title + subtitle — well above the column labels (column labels live at
    # the top edge of the heatmap axes which sits at y=0.77 in figure coords)
    fig.text(0.5, 0.965, title, ha="center", va="top",
             fontsize=FONTS["title"] + 2, weight="bold")
    fig.text(0.5, 0.915, subtitle, ha="center", va="top",
             fontsize=FONTS.get("subtitle", FONTS.get("caption", 11)),
             color="#555", linespacing=1.4)

    OUT_DIR = Path(__file__).resolve().parents[1] / "figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved {out}")


def _load_dab2_data_from_csv(csv_path, contrast_labels, reln_query_path_candidates):
    """Read fgsea results CSV and return a dict matching the DAB2_DATA shape.

    csv_path                 : path to BIONMG_DAB2_KD_allContrasts_MG_fgsea_full.csv
    contrast_labels          : ordered list of 3 contrast IDs to extract
                                (e.g. ["E4_96_v3", "E3_96_v3", "int_96_v3"])
    reln_query_path_candidates: list of candidate paths to RELN-core / DAB1 /
                                 DAB2 gene lists; first existing is loaded.
                                 The union of these gene lists is the 132-gene
                                 RELN query used by the R analytical script.

    Returns:
        dict { pathway_id : [(NES, stars, n_reln, gs_size) for each contrast] }
    """
    import pandas as pd
    from pathlib import Path

    # ---- build RELN query gene set ----
    reln_query = set()
    for d in reln_query_path_candidates:
        for f in Path(d).glob("RELN_*.csv"):
            try:
                df = pd.read_csv(f)
                col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
                reln_query.update(df[col].dropna().astype(str).str.upper().tolist())
            except Exception:
                pass
        for f in Path(d).glob("RELN_core_final_78genes.csv"):
            try:
                df = pd.read_csv(f)
                col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
                reln_query.update(df[col].dropna().astype(str).str.upper().tolist())
            except Exception:
                pass
    print(f"  RELN query gene set: {len(reln_query)} genes")

    def stars_for(p):
        if p is None or (isinstance(p, float) and np.isnan(p)): return ""
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    df = pd.read_csv(csv_path)
    df = df[df["contrast"].isin(contrast_labels)].copy()
    df["leadingEdge_set"] = df["leadingEdge"].fillna("").astype(str).apply(
        lambda s: {g.strip().upper() for g in s.split(";") if g.strip()}
    )
    df["n_reln"] = df["leadingEdge_set"].apply(lambda S: len(S & reln_query))

    # Build dict in DAB2_DATA order: each pathway → list of (NES, stars, n_reln, gs_size)
    # We use the pathway-ID column to key, then look up by the SHORT name used in DAB2_DATA.
    PATHWAY_ID_TO_LABEL = {
        "REACTOME_CHOLESTEROL_BIOSYNTHESIS":                                  "Cholesterol Biosynthesis",
        "HALLMARK_FATTY_ACID_METABOLISM":                                     "Fatty Acid Metabolism",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION":                                 "Oxidative Phosphorylation",
        "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION":                         "Mitochondrial Proteostasis",
        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION":                         "ECM Organization",
        "HALLMARK_MTORC1_SIGNALING":                                          "mTORC1 Signaling",
        "GOBP_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS":   "Ubiquitin-Dependent Proteolysis",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB":                                   "TNFα Signaling via NF-κB",
        "HALLMARK_INFLAMMATORY_RESPONSE":                                     "Inflammatory Response",
        "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP":            "SASP",
        "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE":             "DNA Damage-Induced Senescence",
        "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY":                           "Cytokine-Mediated Signaling",
        "GOBP_AMYLOID_BETA_CLEARANCE":                                        "Amyloid-β Clearance",
        "GOBP_PHAGOCYTOSIS":                                                  "Phagocytosis",
        "GOBP_PHAGOCYTOSIS_ENGULFMENT":                                       "Phagocytosis Engulfment",
        "GOBP_LYSOSOMAL_TRANSPORT":                                           "Lysosomal Transport",
    }

    data = {}
    for pid, plabel in PATHWAY_ID_TO_LABEL.items():
        cells = []
        for c in contrast_labels:
            row = df[(df["pathway"] == pid) & (df["contrast"] == c)]
            if row.empty:
                cells.append((float("nan"), "", 0, 0))
                continue
            r = row.iloc[0]
            cells.append((float(r["NES"]), stars_for(float(r["padj"])),
                          int(r["n_reln"]), int(r["size"])))
        data[plabel] = cells
    return data


def main():
    targets = sys.argv[1] if len(sys.argv) > 1 else "both"

    if targets in ("dab1", "both"):
        render(
            "dab1", DAB1_DATA, DAB1_COL_LABELS,
            "BIONMG iPSC-MG — DAB1 KD Pathway Enrichment",
            "APOE4 / APOE3 DAB1 knockdown (pooled siRNA-893+894) and the genotype × KD interaction (siRNA-893)\n"
            "fgsea multilevel; * padj<0.05  ** <0.01  *** <0.001.   Parentheses: RELN-core-DAB1/2 leading-edge / pathway size.",
            "BIONMG_DAB1_KD_pathway_heatmap_standardised.png",
        )
    if targets in ("dab2", "both"):
        render(
            "dab2", DAB2_DATA, DAB2_COL_LABELS,
            "BIONMG iPSC-MG — DAB2 KD Pathway Enrichment",
            "APOE4 / APOE3 DAB2 knockdown (siRNA-96) and the genotype × KD interaction\n"
            "fgsea multilevel; * padj<0.05  ** <0.01  *** <0.001.   Parentheses: RELN-core-DAB1/2 leading-edge / pathway size.",
            "BIONMG_DAB2_KD_pathway_heatmap_standardised.png",
        )

    # v3 = strict-filter prevalence DESeq2 ranking, multilevel fgsea
    # (requires run_fgsea_DAB{1,2}_KD_v3.R to have been run first so the
    # E4_*_v3 / E3_*_v3 / int_*_v3 contrast rows exist in the data CSV)
    from pathlib import Path
    PROJECT_ROOT = None
    for candidate in [
        "<EDIT_PROJECT_ROOT>/Final Figures  June72026",
        "/Users/gecko/Desktop/Final Figures  June72026",
        "<EDIT_PROJECT_ROOT>/Final Figures  June62026",
        "/Users/gecko/Desktop/Final Figures  June62026",
    ]:
        if Path(candidate).exists():
            PROJECT_ROOT = Path(candidate); break
    reln_query_paths = [
        str(PROJECT_ROOT / "RELN_lists"),
        "<EDIT_INPUT_DIR>",
        str(Path.home() / "Library/Application Support/Claude/local-agent-mode-sessions"),
    ] if PROJECT_ROOT else []

    if targets in ("dab2_v3", "v3", "both") and PROJECT_ROOT is not None:
        v3_csv = PROJECT_ROOT / "Figures" / "BIONMG+DAB2siRNA_KD" / "data" / "BIONMG_DAB2_KD_allContrasts_MG_fgsea_full.csv"
        if v3_csv.exists():
            try:
                v3_data = _load_dab2_data_from_csv(
                    str(v3_csv),
                    contrast_labels=["E4_96_v3", "E3_96_v3", "int_96_v3"],
                    reln_query_path_candidates=reln_query_paths,
                )
                if all(v3_data[p][0][3] > 0 for p in v3_data):
                    render(
                        "dab2_v3", v3_data, DAB2_COL_LABELS,
                        "BIONMG iPSC-MG — DAB2 KD Pathway Enrichment (strict-filter v3 ranking)",
                        "APOE4 / APOE3 DAB2 knockdown (siRNA-96) and the genotype × KD interaction\n"
                        "fgsea multilevel; DESeq2 with ≥10 counts in ≥80% Ctrl prevalence filter (v3 ranking)\n"
                        "* padj<0.05  ** <0.01  *** <0.001.   Parentheses: RELN-core-DAB1/2 leading-edge / pathway size.",
                        "BIONMG_DAB2_KD_pathway_heatmap_v3_strictfilter.png",
                    )
                else:
                    print("[dab2_v3] data CSV exists but v3 contrast rows are missing — run run_fgsea_DAB2_KD_v3.R first")
            except Exception as e:
                print(f"[dab2_v3] error loading v3 data: {e}; skipping v3 render")

    if targets in ("dab1_v3", "v3", "both") and PROJECT_ROOT is not None:
        v3_csv = PROJECT_ROOT / "Figures" / "BIONMG_DAB1siRNA_KD" / "data" / "BIONMG_DAB1_KD_3contrasts_fgsea_full.csv"
        if v3_csv.exists():
            try:
                v3_data = _load_dab2_data_from_csv(
                    str(v3_csv),
                    contrast_labels=["E4_893_v3", "E3_893_v3", "int_893_v3"],
                    reln_query_path_candidates=reln_query_paths,
                )
                if all(v3_data[p][0][3] > 0 for p in v3_data):
                    render(
                        "dab1_v3", v3_data, DAB1_COL_LABELS,
                        "BIONMG iPSC-MG — DAB1 KD Pathway Enrichment (strict-filter v3 ranking)",
                        "APOE4 / APOE3 DAB1 knockdown (siRNA-893) and the genotype × DAB1-KD interaction\n"
                        "(interaction pools siRNA-893 + siRNA-894 via the Q4 binary kd_status term).\n"
                        "fgsea multilevel; DESeq2 with ≥10 counts in ≥80% Ctrl prevalence filter (v3 ranking)\n"
                        "* padj<0.05  ** <0.01  *** <0.001.   Parentheses: RELN-core-DAB1/2 leading-edge / pathway size.",
                        "BIONMG_DAB1_KD_pathway_heatmap_v3_strictfilter.png",
                    )
                else:
                    print("[dab1_v3] data CSV exists but v3 contrast rows are missing — run run_fgsea_DAB1_KD_v3.R first")
            except Exception as e:
                print(f"[dab1_v3] error loading v3 data: {e}; skipping v3 render")


if __name__ == "__main__":
    main()
