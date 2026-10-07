#!/usr/bin/env Rscript
# ============================================================================
# UNBIASED fgseaMultilevel over the FULL Hallmark + KEGG + Reactome universe,
# strict-filter v3, all 7 BIONMG contrasts. RUN THIS IN CURSOR / RStudio.
#
# This replaces the earlier gseapy.prerank (1,000-permutation) unbiased scan with
# a genuine fgseaMultilevel run (adaptive multilevel p-values, nPermSimple=10000),
# matching the rest of the project's canonical fgsea convention. Output feeds
# render_unbiased_HKR_v3_MULTILEVEL_heatmaps.py to draw the heatmaps.
#
# fgsea convention (identical to run_fgsea_unbiased_DAB1_DAB2_v3.R):
#   fgseaMultilevel(minSize=15, maxSize=500, eps=0, nPermSimple=10000,
#                   scoreType="std", seed=42)
#
# Gene sets: MSigDB Hallmark (H) + KEGG legacy + Reactome, via msigdbr
#   (GO:BP intentionally excluded, as in the v3 pipeline). This reproduces the
#   same H+KEGG+Reactome universe the gseapy scan used, but with multilevel stats.
#   NOTE: msigdbr's KEGG_LEGACY should match the c2.cp.kegg_legacy GMT closely;
#   if you need byte-exact gene sets, swap the msigdbr block for read.gmt() on
#   your h.all / c2.cp.reactome / c2.cp.kegg_legacy v2024.1 .gmt files.
#
# Ranked lists: the strict-v3 Wald .rnk files (one per contrast).
# ============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (pkg in c("fgsea", "msigdbr")) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
  for (pkg in c("dplyr", "data.table")) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(dplyr); library(fgsea); library(msigdbr); library(data.table)
})
set.seed(42)

# ---- PATHS: edit these two if your drive layout differs -------------------
RANK_DIR <- "/Volumes/VERBATIM HD/FINAL FIGURES & SCRIPTS_JUNE2026/RELN_v2_redo_Wald/_archived_v3_strict80pct_filter/ranked_Wald_v3_E3E4"
OUT_CSV  <- "/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal/Revised Figures_August2026/Unbiased_HKR_v3_MULTILEVEL/data/MSigDB_HKR_v3_MULTILEVEL_all7contrasts.csv"
# ---------------------------------------------------------------------------
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

# contrast_id -> .rnk filename
CONTRASTS <- c(
  BIONMG_E4_vs_E3_Ctrl               = "BIONMG_E4_vs_E3_Ctrl.rnk",
  BIONMG_DAB1_KD_vs_Ctrl_E3          = "BIONMG_DAB1_KD_vs_Ctrl_E3.rnk",
  BIONMG_DAB1_KD_vs_Ctrl_E4          = "BIONMG_DAB1_KD_vs_Ctrl_E4.rnk",
  BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4 = "BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4.rnk",
  BIONMG_DAB2_KD_vs_Ctrl_E3          = "BIONMG_DAB2_KD_vs_Ctrl_E3.rnk",
  BIONMG_DAB2_KD_vs_Ctrl_E4          = "BIONMG_DAB2_KD_vs_Ctrl_E4.rnk",
  BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4 = "BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4.rnk"
)

# ---- gene sets: Hallmark + KEGG legacy + Reactome -------------------------
cat("[1] loading Hallmark + KEGG + Reactome via msigdbr ...\n")
h  <- msigdbr(species = "Homo sapiens", category = "H")
c2 <- msigdbr(species = "Homo sapiens", category = "C2")
sub <- if ("gs_subcat" %in% names(c2)) c2$gs_subcat else c2$gs_subcollection
reac <- c2[sub == "CP:REACTOME", ]
kegg <- c2[grepl("^CP:KEGG", sub) & !grepl("MEDICUS", sub), ]   # KEGG legacy
gs_df <- bind_rows(
  h[, c("gs_name", "gene_symbol")],
  reac[, c("gs_name", "gene_symbol")],
  kegg[, c("gs_name", "gene_symbol")]
)
gs_list <- split(toupper(gs_df$gene_symbol), gs_df$gs_name)
cat(sprintf("    H=%d  KEGG=%d  Reactome=%d  total=%d sets\n",
            length(unique(h$gs_name)), length(unique(kegg$gs_name)),
            length(unique(reac$gs_name)), length(gs_list)))

read_rnk <- function(path) {
  d <- fread(path, header = FALSE, col.names = c("gene", "stat"))
  d <- d[!is.na(gene) & !is.na(stat)]
  d[, gene := toupper(as.character(gene))]
  d <- d[order(-abs(stat))][, .SD[1], by = gene][order(-stat)]
  setNames(d$stat, d$gene)
}

collection_of <- function(p) ifelse(startsWith(p, "HALLMARK_"), "Hallmark",
                             ifelse(startsWith(p, "KEGG_"), "KEGG",
                             ifelse(startsWith(p, "REACTOME_"), "Reactome", "Other")))

run_one <- function(rnk_path, cid) {
  set.seed(42)
  stats <- read_rnk(rnk_path)
  res <- fgseaMultilevel(pathways = gs_list, stats = stats,
                         minSize = 15, maxSize = 500, eps = 0,
                         nPermSimple = 10000, scoreType = "std",
                         BPPARAM = BiocParallel::SerialParam(), nproc = 1)
  setDT(res)
  res[, contrast_id := cid]
  res[, collection  := collection_of(pathway)]
  res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  res[, padj := p.adjust(pval, method = "BH")]           # BH per contrast
  res[, .(contrast_id, pathway, collection, ES, NES, pval, padj, size, leadingEdge)]
}

all_res <- rbindlist(lapply(names(CONTRASTS), function(cid) {
  p <- file.path(RANK_DIR, CONTRASTS[[cid]])
  cat(sprintf("[2] fgseaMultilevel: %s\n", cid))
  run_one(p, cid)
}), use.names = TRUE)

fwrite(all_res, OUT_CSV)
cat(sprintf("\nDONE. wrote %s (%d rows across %d contrasts)\n",
            OUT_CSV, nrow(all_res), length(CONTRASTS)))
cat("Next: python3 render_unbiased_HKR_v3_MULTILEVEL_heatmaps.py\n")
