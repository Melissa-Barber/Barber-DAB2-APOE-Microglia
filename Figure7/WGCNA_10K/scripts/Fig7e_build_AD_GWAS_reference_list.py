#!/usr/bin/env python3
"""
Build the AD-GWAS reference panel directly from the Bellenguez et al. 2022
Nat Genet supplementary workbook (10.1038/s41588-022-01024-z).

All 75 risk loci:
  * 35 known-locus genes            -> Supplementary Table 5 ("Gene" col, rows
                                       whose "Known locus" label is NOT "New")
  * 55 new-locus prioritised genes  -> Supplementary Table 20 (Gene Prioritization
                                       Tier == "Tier 1" or "Tier 2")
  * APOE                            -> added; the APOE region is the established
                                       dominant signal and is NOT among the 75 loci.

Writes AD_GWAS_reference_list_Bellenguez2022.csv (gene_symbol, locus_status, tier,
source). Run this BEFORE render_AD_GWAS_enrichment_Bellenguez2022.py.
"""
from pathlib import Path
import openpyxl
import pandas as pd

# self-contained: the supplementary xlsx is kept beside this script's bundle
XLSX = Path(__file__).resolve().parent / "Bellenguez2022_supplementary_tables.xlsx"
if not XLSX.exists():   # fallback to the WGCNA bundle copy
    XLSX = Path("<EDIT_PROJECT_ROOT>/Final Figures /"
                "Bellengeuz et al 2022_Table 5 all 91 loci AD risk microglia/"
                "Bellenguez2022_supplementary_tables.xlsx")
OUT = Path("<EDIT_PROJECT_ROOT>/Final Figures /WGCNA_10K/Module Gene lists")

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)

# ---- Supplementary Table 5: known-locus genes ----
r5 = list(wb["Supplementary Table 5"].iter_rows(values_only=True))
known = set()
for r in r5[4:]:                       # header rows 0-3
    if r[0] in (None, "") or r[3] in (None, ""):
        continue
    locus_label = str(r[4]).strip() if r[4] else ""
    if locus_label.lower() != "new":   # known loci are labelled by their locus name
        known.add(str(r[3]).strip())

# ---- Supplementary Table 20: new-locus tier 1 / tier 2 prioritised genes ----
r20 = list(wb["Supplementary Table 20"].iter_rows(values_only=True))
t1, t2 = set(), set()
for r in r20[3:]:                      # header rows 0-2
    if r[0] in (None, ""):
        continue
    tier = str(r[8]).strip() if r[8] else ""
    if tier == "Tier 1":
        t1.add(str(r[0]).strip())
    elif tier == "Tier 2":
        t2.add(str(r[0]).strip())

rows = ([(g, "Known", "known-lead") for g in sorted(known)] +
        [(g, "New", "Tier 1") for g in sorted(t1)] +
        [(g, "New", "Tier 2") for g in sorted(t2)] +
        [("APOE", "APOE-region", "added (not in the 75 loci; established signal)")])
panel = pd.DataFrame(rows, columns=["gene_symbol", "locus_status", "tier"]).drop_duplicates("gene_symbol")
panel = panel[~panel.gene_symbol.str.contains(" ", na=False)]      # drop "IGH gene cluster" etc.
panel["source"] = "Bellenguez 2022 NatGen 10.1038/s41588-022-01024-z | SuppT5 known-lead + SuppT20 Tier1/2 (+APOE)"
panel.to_csv(OUT / "AD_GWAS_reference_list_Bellenguez2022.csv", index=False)
print(f"panel: {len(panel)} genes "
      f"(known-lead {(panel.tier=='known-lead').sum()}, Tier1 {(panel.tier=='Tier 1').sum()}, "
      f"Tier2 {(panel.tier=='Tier 2').sum()}, +APOE) -> AD_GWAS_reference_list_Bellenguez2022.csv")
