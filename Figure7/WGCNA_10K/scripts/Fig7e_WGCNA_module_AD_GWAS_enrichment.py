#!/usr/bin/env python3
"""
AD-GWAS risk-gene enrichment per WGCNA-10K module — AUTHORITATIVE reference list
from Bellenguez et al. 2022 Nat Genet (10.1038/s41588-022-01024-z), all 75 loci:
  - 35 known-locus genes            (Supplementary Table 5, "Gene" col, known loci)
  - 55 new-locus prioritised genes  (Supplementary Table 20, Tier 1 + Tier 2)
  - APOE                            (added; the APOE region is the established
                                     dominant signal and is NOT among the 75 loci)
Self-contained: reads the saved reference list + the WGCNA master gene table,
runs a one-sided hypergeometric test per module (BH-FDR), writes the stats CSV
and renders the bar figure.
"""
from pathlib import Path
import sys
import pandas as pd
from scipy.stats import hypergeom
import matplotlib.pyplot as plt

ROOT = Path("<EDIT_PROJECT_ROOT>/Final Figures /WGCNA_10K")
GLST = ROOT / "Module Gene lists"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import apply_style, CM, DPI
import module_labels_10K as L
apply_style()

# ---- reference panel + WGCNA modules ----
panel = pd.read_csv(GLST / "AD_GWAS_reference_list_Bellenguez2022.csv")
ref = set(panel["gene_symbol"].astype(str).str.upper())
g = pd.read_csv(ROOT / "tables" / "WGCNA_10K_master_gene_table_p12.csv")
g["sym"] = g["gene_symbol"].astype(str).str.upper()
present = g[g["sym"].isin(ref)]
N, K = len(g), len(present)

STATE = {"brown": "CRM", "green": "no state", "magenta": "HM", "blue": "RM/HLA",
         "greenyellow": "DAM", "salmon": "HM", "purple": "DAM/LDAM"}
rows = []
for m, n in g["module"].value_counts().items():
    k = int((present["module"] == m).sum())
    fe = (k / n) / (K / N) if k else 0.0
    p = hypergeom.sf(k - 1, N, K, n) if k else 1.0
    rows.append(dict(module=m, mg_state=STATE.get(m, "unassigned" if m == "grey" else "n.s."),
                     n_module=n, k_AD=k, expected=round(n * K / N, 2),
                     fold_enrichment=round(fe, 2), p_hyper=p))
res = pd.DataFrame(rows).sort_values("p_hyper").reset_index(drop=True)
res["q_BH"] = (res["p_hyper"] * len(res) / (res.index + 1)).clip(upper=1.0)
res["q_BH"] = res["q_BH"][::-1].cummin()[::-1]
res["sig"] = res["q_BH"].map(lambda q: "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns")
res.to_csv(GLST / "AD_GWAS_module_enrichment_Bellenguez2022.csv", index=False)

# ---- figure ----
COL = dict(L.MODULE_COLORS_10K)
COL.update({"black": "#333333", "grey": "#9E9E9E", "turquoise": "#1FB6B6",
            "yellow": "#D4B100", "red": "#D62728", "pink": "#E377C2",
            "tan": "#B59A6B", "cyan": "#17A2B8"})
res = res.sort_values(["p_hyper", "fold_enrichment"], ascending=[True, False]).reset_index(drop=True)
n = len(res)
T_MARG, B_MARG = 1.4, 2.6
fh = T_MARG + 0.62 * n + B_MARG
fig = plt.figure(figsize=(20.0 * CM, fh * CM))
ax = fig.add_axes([0.28, B_MARG / fh, 0.56, 1 - (T_MARG + B_MARG) / fh])
fe_max = res["fold_enrichment"].max()
for i, r in res.iterrows():
    y = n - 1 - i
    c = COL.get(r["module"], "#888")
    ax.barh(y, r["fold_enrichment"], height=0.66, color=c, edgecolor="#333", linewidth=0.5, zorder=3)
    star = "" if r["sig"] == "ns" else r["sig"]
    qtxt = f"{r['q_BH']:.0e}".replace("e-0", "e-") if r["q_BH"] < 0.01 else f"{r['q_BH']:.2f}"
    ax.text(r["fold_enrichment"] + fe_max * 0.015, y, f"  {r['k_AD']}/{r['n_module']}   q={qtxt} {star}",
            va="center", ha="left", fontsize=8, color="#222", fontweight="bold" if star else "normal")
    state = "" if r["mg_state"] in ("n.s.", "unassigned") else f" — {r['mg_state']}"
    ax.text(-0.02, y, f"{r['module']}{state}", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=9, fontweight="bold", color=c, clip_on=False)
ax.axvline(1.0, color="#888", lw=0.9, ls="--", zorder=1)
ax.text(1.0, n - 0.35, "no enrichment", fontsize=6.5, color="#888", ha="center", va="bottom")
ax.set_xlim(0, fe_max * 1.42); ax.set_ylim(-0.6, n - 0.4); ax.set_yticks([])
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
ax.set_xlabel("Fold enrichment (observed / expected AD-GWAS genes)", fontsize=9, labelpad=5)
ax.tick_params(axis="x", labelsize=8); ax.grid(axis="x", color="#e3e3e3", lw=0.5, zorder=0)
fig.suptitle("AD-GWAS risk-gene enrichment per WGCNA-10K module  (Bellenguez 2022, all 75 loci)",
             x=0.012, y=0.985, ha="left", fontsize=11.5, fontweight="bold")
fig.text(0.012, 0.008,
         f"One-sided hypergeometric test; background = the {N:,}-gene WGCNA network; "
         f"K={K} reference genes present ({len(ref)-K} of {len(ref)} not in the top-10K variable set). "
         "Bars = fold enrichment; k/n = AD-GWAS genes / module size; q = BH-FDR; *q<0.05 **q<0.01. "
         "Reference = Bellenguez 2022 (EADB) all 75 loci: 35 known-locus genes (Suppl. Table 5) + 55 "
         "new-locus Tier 1/2 prioritised genes (Suppl. Table 20) + APOE (APOE region not among the 75 loci).",
         ha="left", va="bottom", fontsize=6, color="#555", wrap=True)
out = ROOT / "figures" / "Fig_10K_AD_GWAS_module_enrichment_Bellenguez2022.png"
fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"saved {out.name}  | N={N} K={K} | sig:",
      ", ".join(f"{r.module}(q={r.q_BH:.3f})" for _, r in res[res.sig != 'ns'].iterrows()))
