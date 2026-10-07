#!/usr/bin/env python3
"""
GROUP-AWARE pipeline — step 3  (RANK BUILDER ONLY — no fgsea here).

Builds the Wald-stat ranked .rnk files from the group-aware DESeq2 CSVs, mapping
ensembl -> hgnc_symbol via the existing DEG mapping file. The .rnk file names are
IDENTICAL to the strict-filter ones (BIONMG_E4_vs_E3_Ctrl.rnk, etc.), so the
fgseaMultilevel R script picks them up unchanged — you only repoint RANK_ARCH at
this folder.

IMPORTANT: fgsea is NOT run in Python. This script deliberately stops at the .rnk
files. All enrichment scoring is done by run_fgsea_MULTILEVEL_v3_ALL.R
(fgseaMultilevel, nPermSimple=10000, eps=0, scoreType="std", seed=42). Do not add a
gseapy / prerank step here — keeping the permutation engine in R is the whole reason
the results are stable.
"""
from pathlib import Path
import pandas as pd
warnings = __import__("warnings"); warnings.filterwarnings("ignore")

# ---- PATHS: must match steps 01/02 ----
ROOT = Path("/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal")
GA   = ROOT / "Revised Figures_August2026" / "GroupAware_filter_pipeline"
V3   = GA / "DESeq2_groupaware"                       # step-02 output CSVs
RANK = GA / "ranked_Wald_groupaware_E3E4"; RANK.mkdir(parents=True, exist_ok=True)

# ensembl -> symbol map (reuse the same DEG mapping used by the strict pipeline)
deg = pd.read_csv(ROOT / "DEG" / "Q1_Genotype_Ctrl" / "E2_vs_E3_Ctrl.csv")
sym_map = dict(zip(deg["ensembl_id"], deg["hgnc_symbol"]))

plan = {
    "BIONMG_E4_vs_E3_Ctrl":                 "Q1_E4_vs_E3_Ctrl.csv",
    "BIONMG_DAB1_KD_vs_Ctrl_E3":            "Q2_DAB1_KD_vs_Ctrl_E3.csv",
    "BIONMG_DAB1_KD_vs_Ctrl_E4":            "Q2_DAB1_KD_vs_Ctrl_E4.csv",
    "BIONMG_DAB2_KD_vs_Ctrl_E3":            "Q3_DAB2_KD_vs_Ctrl_E3.csv",
    "BIONMG_DAB2_KD_vs_Ctrl_E4":            "Q3_DAB2_KD_vs_Ctrl_E4.csv",
    "BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4":   "Q3d_DAB1_KD_vs_Ctrl_pooled_E3E4.csv",
    "BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4":   "Q3d_DAB2_KD_vs_Ctrl_pooled_E3E4.csv",
}

for cid, fn in plan.items():
    d = pd.read_csv(V3 / fn).dropna(subset=["stat"])
    d["gene"] = d["ensembl_id"].map(sym_map)
    d = d.dropna(subset=["gene"])[["gene", "stat"]]
    d["abs"] = d.stat.abs()
    d = (d.sort_values("abs", ascending=False)
           .drop_duplicates("gene")
           .drop(columns="abs")
           .sort_values("stat", ascending=False))
    d.to_csv(RANK / f"{cid}.rnk", sep="\t", header=False, index=False)
    print(f"{cid}.rnk  ({len(d)} genes)")

print(f"\nwrote {len(plan)} .rnk files to {RANK}")
print("next: set RANK_ARCH in run_fgsea_MULTILEVEL_v3_ALL.R to this folder, then run it.")
