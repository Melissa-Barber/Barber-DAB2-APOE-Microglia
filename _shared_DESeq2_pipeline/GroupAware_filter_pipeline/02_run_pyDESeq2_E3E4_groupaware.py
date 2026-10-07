#!/usr/bin/env python3
"""
GROUP-AWARE pipeline — step 2  (identical modelling to the strict-filter step 2,
only the input counts file and the output directory change).

Runs pyDESeq2 with a Batch_grp covariate for the priority E3/E4 contrasts:
  Q1_E4_vs_E3_Ctrl
  Q2_DAB1_KD_vs_Ctrl_{E3,E4}
  Q3_DAB2_KD_vs_Ctrl_{E3,E4}
  Q3d_{DAB1,DAB2}_KD_vs_Ctrl_pooled_E3E4 (genotype-adjusted)
Each result CSV has columns: ensembl_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.

NOTHING about the statistics changes vs the strict run — same design, same Wald
test, same contrasts — so any difference in the fgsea result is attributable ONLY
to the gene-detection filter (group-aware vs Ctrl-only). That is the whole point of
keeping this as a parallel pipeline.
"""
from pathlib import Path
import pandas as pd, warnings
warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

# ---- PATHS: must match step 01 ----
ROOT = Path("/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal")
OUT  = ROOT / "Revised Figures_August2026" / "GroupAware_filter_pipeline" / "DESeq2_groupaware"

counts = pd.read_csv(OUT / "counts_filtered_groupaware80pct.csv", index_col=0)
meta   = pd.read_csv(OUT / "sample_metadata.csv", index_col=0)
meta.index = meta.index.astype(str); counts.columns = counts.columns.astype(str)
counts_T = counts.T.astype(int)
bc = meta["Batch"].value_counts(); keep = bc[bc >= 3].index.tolist()
meta["Batch_grp"] = meta["Batch"].where(meta["Batch"].isin(keep), "Other")
inf = DefaultInference(n_cpus=4)


def run(mask, design, contrast, out_name):
    sm = meta.loc[mask].copy(); sc = counts_T.loc[sm.index]
    df = [f for f in design if sm[f].nunique() >= 2]
    dds = DeseqDataSet(counts=sc, metadata=sm, design_factors=df, inference=inf, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=list(contrast), inference=inf, quiet=True); st.summary()
    res = st.results_df.copy(); res.index.name = "ensembl_id"
    res.to_csv(OUT / f"{out_name}.csv")
    print(f"{out_name}: n={len(sm)} sig padj<0.05 = {(res.padj < 0.05).sum()}")


run((meta.siRNA == "Ctrl") & meta.Genotype.isin(["E3", "E4"]),
    ["Batch_grp", "Genotype"], ["Genotype", "E4", "E3"], "Q1_E4_vs_E3_Ctrl")
for g in ["E3", "E4"]:
    run((meta.Genotype == g) & meta.KD_class.isin(["Ctrl", "DAB1"]),
        ["Batch_grp", "KD_class"], ["KD_class", "DAB1", "Ctrl"], f"Q2_DAB1_KD_vs_Ctrl_{g}")
for g in ["E3", "E4"]:
    run((meta.Genotype == g) & meta.KD_class.isin(["Ctrl", "DAB2"]),
        ["Batch_grp", "KD_class"], ["KD_class", "DAB2", "Ctrl"], f"Q3_DAB2_KD_vs_Ctrl_{g}")
for kd in ["DAB1", "DAB2"]:
    run(meta.Genotype.isin(["E3", "E4"]) & meta.KD_class.isin(["Ctrl", kd]),
        ["Batch_grp", "Genotype", "KD_class"], ["KD_class", kd, "Ctrl"],
        f"Q3d_{kd}_KD_vs_Ctrl_pooled_E3E4")
