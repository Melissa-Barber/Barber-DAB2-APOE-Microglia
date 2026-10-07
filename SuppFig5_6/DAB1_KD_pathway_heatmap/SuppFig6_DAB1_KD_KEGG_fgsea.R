#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# DAB1 siRNA KD — KEGG (legacy) fgsea pass.
#
# Direct sibling of the existing DAB2 KEGG file
#   Figures/BIONMG+DAB2siRNA_KD/data/BIONMG_DAB2_KD_KEGG_fgsea_full.csv
# and of the Ctrl file BIONMG_E4vsE3_KEGG_fgsea_multilevel.csv.
# The DAB1 KD unbiased fgsea was originally run on Hallmark + Reactome + GO:BP
# (companion file geneset_membership_H_GOBP_Reactome.csv) with KEGG excluded;
# this fills that gap for DAB1, matching the KEGG convention used elsewhere.
#
# Recipe (identical to fill_gaps_fgsea_multilevel.R section (3)):
#   fgsea::fgsea(pathways = KEGG_legacy_GMT, stats,
#                minSize = 15, maxSize = 500, eps = 0,
#                nPermSimple = 10000, scoreType = "std", seed = 42)
#   BH per contrast; read_rnk dedups by max |stat|.
#   Gene sets: c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt (186 sets), uppercased.
#
# Contrasts mirror the DAB2 KEGG file (E3/E4/pool per siRNA + interaction):
#   E3_893, E4_893, pool_893, E3_894, E4_894, pool_894, int_893
#
# Output schema matches the DAB2 KEGG file exactly:
#   pathway, pval, padj, log2err, ES, NES, size, leadingEdge, contrast
#
# How to run: open in Cursor / RStudio and Source, or
#   Rscript "run_fgsea_DAB1_KD_KEGG.R"
# Tested with fgsea 1.30+, R 4.3+.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("fgsea", quietly = TRUE))       BiocManager::install("fgsea", update = FALSE, ask = FALSE)
  library(fgsea); library(data.table)
})
set.seed(42)

# ---- PATHS (edit only if your drive layout differs) -----------------------
PROJECT <- "/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal"
RANK    <- file.path(PROJECT, "ranked_Wald")
GMT     <- file.path(PROJECT, "Figures", "MsigDB_GeneMatrixTransposition",
                     "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")
OUT     <- file.path(PROJECT, "Revised Figures September 2026", "data", "SuppFig5_6",
                     "DAB1_KD_pathway_heatmap", "BIONMG_DAB1_KD_KEGG_fgsea_full.csv")
# ---------------------------------------------------------------------------
stopifnot(file.exists(GMT))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

# contrast label -> .rnk file (mirrors the DAB2 KEGG contrast set)
MAP <- c(
  E3_893   = "BIONMG_893_vs_Ctrl_E3.rnk",
  E4_893   = "BIONMG_893_vs_Ctrl_E4.rnk",
  pool_893 = "BIONMG_893_vs_Ctrl_pooled.rnk",
  E3_894   = "BIONMG_894_vs_Ctrl_E3.rnk",
  E4_894   = "BIONMG_894_vs_Ctrl_E4.rnk",
  pool_894 = "BIONMG_894_vs_Ctrl_pooled.rnk",
  int_893  = "BIONMG_DAB1_response_E4_vs_E3.rnk"
)

read_rnk <- function(f) {
  d <- fread(f, header = FALSE, col.names = c("g", "s"))
  d <- d[!is.na(g) & !is.na(s)]
  d[, abs_s := abs(s)]
  d <- d[order(-abs_s)][, .SD[1], by = g][order(-s)]
  setNames(d$s, toupper(d$g))
}

kegg_gmt <- gmtPathways(GMT)
kegg_gmt <- lapply(kegg_gmt, function(g) unique(toupper(g)))
cat(sprintf("KEGG legacy gene sets: %d\n", length(kegg_gmt)))

run_one <- function(stats, contrast) {
  r <- fgsea(pathways = kegg_gmt, stats = stats,
             minSize = 15, maxSize = 500, eps = 0,
             nPermSimple = 10000, scoreType = "std",
             BPPARAM = BiocParallel::SerialParam(), nproc = 1)
  if (nrow(r) == 0) return(data.table())
  setnames(r, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge"))
  r[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  r[, padj := p.adjust(pval, method = "BH")]   # BH per contrast
  r[, contrast := contrast]
  r[, .(pathway, pval, padj, log2err, ES, NES, size, leadingEdge, contrast)]
}

all_res <- rbindlist(lapply(names(MAP), function(cid) {
  f <- file.path(RANK, MAP[[cid]])
  stopifnot(file.exists(f))
  cat(sprintf("  fgsea KEGG: %-9s <- %s\n", cid, MAP[[cid]]))
  run_one(read_rnk(f), cid)
}), use.names = TRUE)

fwrite(all_res, OUT)
cat(sprintf("\nDONE. wrote %s\n  (%d rows, %d KEGG pathways x %d contrasts)\n",
            OUT, nrow(all_res), all_res[, uniqueN(pathway)], length(MAP)))
# quick peek: complement / inflammation KEGG sets in the pooled DAB1-KD columns
peek <- all_res[grepl("COMPLEMENT|INFLAMMAT|LYSOSOME|PHAGO|APOPTOSIS", pathway) &
                contrast %in% c("pool_893", "E4_893", "E3_893")]
if (nrow(peek)) print(peek[order(contrast, padj), .(contrast, pathway, NES, padj, size)])
