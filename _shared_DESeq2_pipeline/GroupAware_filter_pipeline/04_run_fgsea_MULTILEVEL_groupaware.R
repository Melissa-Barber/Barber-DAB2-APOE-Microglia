#!/usr/bin/env Rscript
# ============================================================================
# GROUP-AWARE multilevel-fgsea run — RUN IN CURSOR / RStudio (R, not Python).
#
# Identical engine to the strict master (run_fgsea_MULTILEVEL_v3_ALL.R):
#   fgseaMultilevel(nPermSimple=10000, eps=0, scoreType="std", seed=42)
# on the SAME Hallmark + KEGG + Reactome + GO:BP universe. The ONLY difference
# is the input ranks: this reads the group-aware .rnk files (built by
# 03_build_ranks_groupaware.py from the group-aware DESeq2). So any change in the
# result vs the strict run is due to the gene-detection filter and nothing else.
#
# DO NOT change nPermSimple, eps, scoreType, minSize/maxSize, or switch to a
# Python/gseapy engine — keeping these identical to the strict run is what makes
# this a clean filter-only sensitivity comparison.
#
# Runs two sets (the ones the group-aware ranks cover):
#   1. data/MSigDB_HKR_v3_MULTILEVEL_topupdown_7contrasts_GROUPAWARE.csv
#        universe H+KEGG+Reactome+GO:BP (minSize 15, maxSize 500)
#        contrasts: E4-vs-E3; DAB1_KD & DAB2_KD (E3 / E4 / pooled)
#        -> single-column top-15-up/10-down heatmaps
#   3. data/FOCUSED_v3_MULTILEVEL_7contrasts_GROUPAWARE.csv
#        curated 58-pathway MG panel (H+KEGG+Reactome+GO:BP; minSize 5, maxSize 2000)
#        -> focused MG categorized heatmap
#
# SET 2 (3-column per-siRNA interaction) is intentionally NOT run here: the
# per-siRNA (893/894/96/98) ranks are not rebuilt by the group-aware pipeline.
# Add per-siRNA contrasts to steps 02/03 first if you want those figures too.
# ============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in c("fgsea","msigdbr")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
  for (p in c("dplyr","data.table")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(dplyr); library(fgsea); library(msigdbr); library(data.table)
})
set.seed(42)

# ---- PATHS: edit if your drive layout differs -----------------------------
PROJ      <- "/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal"
GA        <- file.path(PROJ, "Revised Figures_August2026/GroupAware_filter_pipeline")
RANK_GA   <- file.path(GA, "ranked_Wald_groupaware_E3E4")   # group-aware ranks (step 03)
OUT_DIR   <- file.path(GA, "data")
# focused panel definition — reuse the bundled one from the strict multilevel folder
PANEL_CSV <- file.path(PROJ, "Revised Figures_August2026/Unbiased_HKR_v3_MULTILEVEL/data/FOCUSED_panel_definition_58.csv")
# ---------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- gene sets (identical to the strict master) ---------------------------
cat("[1] loading MSigDB collections via msigdbr ...\n")
h    <- msigdbr(species="Homo sapiens", category="H")
c2   <- msigdbr(species="Homo sapiens", category="C2")
csub <- if ("gs_subcat" %in% names(c2)) c2$gs_subcat else c2$gs_subcollection
reac <- c2[csub == "CP:REACTOME", ]
kegg <- c2[grepl("^CP:KEGG", csub) & !grepl("MEDICUS", csub), ]           # KEGG legacy
gobp <- msigdbr(species="Homo sapiens", category="C5")
gsub <- if ("gs_subcat" %in% names(gobp)) gobp$gs_subcat else gobp$gs_subcollection
gobp <- gobp[gsub == "GO:BP", ]

mk <- function(df) split(toupper(df$gene_symbol), df$gs_name)
HKR   <- c(mk(h), mk(kegg), mk(reac))
HKRGO <- c(HKR, mk(gobp))
cat(sprintf("    H=%d KEGG=%d Reactome=%d  -> HKR=%d sets  (+GO:BP=%d)\n",
            length(unique(h$gs_name)), length(unique(kegg$gs_name)),
            length(unique(reac$gs_name)), length(HKR), length(unique(gobp$gs_name))))

# ---- helpers (identical to the strict master) -----------------------------
read_rnk <- function(path){
  d <- fread(path, header=FALSE, col.names=c("gene","stat"))
  d <- d[!is.na(gene) & !is.na(stat)][, gene := toupper(as.character(gene))]
  d <- d[order(-abs(stat))][, .SD[1], by=gene][order(-stat)]
  setNames(d$stat, d$gene)
}
coll_of <- function(p) ifelse(startsWith(p,"HALLMARK_"),"Hallmark",
                       ifelse(startsWith(p,"KEGG_"),"KEGG",
                       ifelse(startsWith(p,"REACTOME_"),"Reactome",
                       ifelse(startsWith(p,"GOBP_"),"GO:BP","Other"))))
run_one <- function(rnk_path, cid, sets, minSize, maxSize){
  set.seed(42)
  res <- fgseaMultilevel(pathways=sets, stats=read_rnk(rnk_path),
                         minSize=minSize, maxSize=maxSize, eps=0,
                         nPermSimple=10000, scoreType="std",
                         BPPARAM=BiocParallel::SerialParam(), nproc=1)
  setDT(res)
  res[, contrast_id := cid]
  res[, collection  := coll_of(pathway)]
  res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse=";")]
  res[, padj := p.adjust(pval, method="BH")]
  res[, .(contrast_id, pathway, collection, ES, NES, pval, padj, size, leadingEdge)]
}
run_many <- function(map, rank_dir, sets, minSize, maxSize){
  rbindlist(lapply(names(map), function(cid){
    cat(sprintf("    fgsea: %s\n", cid))
    run_one(file.path(rank_dir, map[[cid]]), cid, sets, minSize, maxSize)
  }), use.names=TRUE)
}

# ---- SET 1: unbiased top-up/down (group-aware ranks) ----------------------
cat("[2] SET 1 — unbiased top-up/down (H+K+R+GO:BP), group-aware ranks ...\n")
A <- c(BIONMG_E4_vs_E3_Ctrl="BIONMG_E4_vs_E3_Ctrl.rnk",
       BIONMG_DAB1_KD_vs_Ctrl_E3="BIONMG_DAB1_KD_vs_Ctrl_E3.rnk",
       BIONMG_DAB1_KD_vs_Ctrl_E4="BIONMG_DAB1_KD_vs_Ctrl_E4.rnk",
       BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4="BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4.rnk",
       BIONMG_DAB2_KD_vs_Ctrl_E3="BIONMG_DAB2_KD_vs_Ctrl_E3.rnk",
       BIONMG_DAB2_KD_vs_Ctrl_E4="BIONMG_DAB2_KD_vs_Ctrl_E4.rnk",
       BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4="BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4.rnk")
fwrite(run_many(A, RANK_GA, HKRGO, 15, 500),
       file.path(OUT_DIR, "MSigDB_HKR_v3_MULTILEVEL_topupdown_7contrasts_GROUPAWARE.csv"))

# ---- SET 3: curated 58-pathway MG focused panel (group-aware ranks) -------
cat("[3] SET 3 — focused MG panel (H+K+R+GO:BP), group-aware ranks ...\n")
panel <- fread(PANEL_CSV)                                  # pathway, category, setSize
focused_sets <- HKRGO[intersect(panel$pathway, names(HKRGO))]
missing <- setdiff(panel$pathway, names(HKRGO))
if (length(missing)) cat(sprintf("    NOTE: %d panel pathways not found in msigdbr: %s\n",
                                  length(missing), paste(missing, collapse=", ")))
focres <- run_many(A, RANK_GA, focused_sets, 5, 2000)
focres <- merge(focres, panel[, .(pathway, category)], by="pathway", all.x=TRUE)
fwrite(focres, file.path(OUT_DIR, "FOCUSED_v3_MULTILEVEL_7contrasts_GROUPAWARE.csv"))

cat("\nDONE. Wrote 2 CSVs to", OUT_DIR,
    "\nSend me the two _GROUPAWARE.csv files and I'll render the heatmaps,",
    "\nor run the render_*.py plotters pointed at these CSVs.\n")
