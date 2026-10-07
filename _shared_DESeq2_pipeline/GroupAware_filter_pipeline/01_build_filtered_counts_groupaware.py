#!/usr/bin/env python3
"""
GROUP-AWARE prevalence filter (step 1) — replaces the Ctrl-only 80% filter.

WHY: the old filter kept genes detected (>=10 counts) in >=80% of *Ctrl-siRNA*
samples only. That silently drops genes that are condition-specific — e.g. a gene
switched ON by a knockdown (off in Ctrl) or expressed in only one genotype — which
is exactly the biology of a KD or genotype contrast. This is why the ε4-vs-ε3
baseline (both arms Ctrl-siRNA, so a gene ON in ε4 / OFF in ε3 is only ~59% of the
pooled Ctrl samples) lost real candidates.

RULE (group-aware, applied uniformly to every contrast):
    keep a gene if it is detected (>=10 counts) in >= FRAC of the samples of
    AT LEAST ONE genotype x siRNA group (groups with < MIN_GROUP_N samples are
    not used as a denominator). This is the standard "expressed in at least one
    condition" filter (cf. edgeR::filterByExpr) — the principled, uniform choice.

Output: counts_filtered_groupaware80pct.csv + sample_metadata.csv, written to a
NEW directory so the strict-filter results are preserved for the sensitivity
comparison. Run this in Cursor, then 02 -> 03 -> the fgseaMultilevel R script.
"""
from pathlib import Path
import pandas as pd, numpy as np, re

# ---- PATHS: edit if your layout differs -----------------------------------
ROOT = Path("/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal")
GC   = ROOT / "GeneCounts"
OUT  = ROOT / "Revised Figures_August2026" / "GroupAware_filter_pipeline" / "DESeq2_groupaware"
# ---------------------------------------------------------------------------
OUT.mkdir(parents=True, exist_ok=True)

FRAC = 0.80          # fraction of a group that must detect the gene
MIN_GROUP_N = 3      # groups smaller than this are not valid denominators
MIN_COUNT = 10       # "detected" = >= this many counts


def read_fc(p):
    d = pd.read_csv(p, sep="\t", comment="#")
    samp = [c for c in d.columns if c not in ["Geneid", "Chr", "Start", "End", "Strand", "Length"]]
    d = d.set_index("Geneid")[samp]
    rename = {c: (re.match(r"^(\w+)_(S\d+)_", c).group(1) if re.match(r"^(\w+)_(S\d+)_", c) else c) for c in samp}
    return d.rename(columns=rename)


c1 = read_fc(GC / "MelissaB.gene.counts.txt")
c2 = read_fc(GC / "gene.counts.txt")
c2.columns = [c.lstrip("S") if re.match(r"^S\d+$", c) else c for c in c2.columns]
joint = c1.join(c2, how="outer").fillna(0).astype(int)
if "Undetermined" in joint.columns:
    joint = joint.drop(columns="Undetermined")

m_main = pd.read_excel(GC / "MetaDataBatch.xlsx"); m_main["Sample ID"] = m_main["Sample ID"].astype(str)
samp_meta = m_main.set_index("Sample ID")
samp_used = sorted(set(joint.columns) & set(samp_meta.index))
counts = joint[samp_used]; meta = samp_meta.loc[samp_used].copy()
meta["KD_class"] = meta["siRNA"].astype(str).map(
    lambda x: "Ctrl" if x == "Ctrl" else ("DAB1" if x in ["893", "894"] else ("DAB2" if x in ["96", "98"] else "other")))

# ---- group-aware filter ----
grp = meta["Genotype"].astype(str) + "|" + meta["siRNA"].astype(str)
gene_pass = pd.Series(False, index=counts.index)
print("group-aware prevalence filter (>=%d counts in >=%.0f%% of >=1 group, groups >=%d samples):"
      % (MIN_COUNT, 100 * FRAC, MIN_GROUP_N))
for gname, cols in grp.groupby(grp).groups.items():
    cols = list(cols); n = len(cols)
    if n < MIN_GROUP_N:
        print(f"   {gname:12} n={n:2}  skipped (<{MIN_GROUP_N})"); continue
    thr = int(np.ceil(FRAC * n))
    passes = (counts[cols] >= MIN_COUNT).sum(axis=1) >= thr
    gene_pass = gene_pass | passes
    print(f"   {gname:12} n={n:2}  thr={thr:2}  genes_detected={int(passes.sum())}")

filt = counts.loc[gene_pass]
filt.to_csv(OUT / "counts_filtered_groupaware80pct.csv")
meta.to_csv(OUT / "sample_metadata.csv")
print(f"\nsaved counts {filt.shape} + metadata {meta.shape} to {OUT}")
print(f"genes passing group-aware filter = {len(filt)}   "
      f"(compare with the strict Ctrl-only filter = 14,642)")
