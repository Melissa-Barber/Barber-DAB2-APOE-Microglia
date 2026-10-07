#!/usr/bin/env Rscript
# ============================================================================
# MASTER multilevel-fgsea run for ALL v3 heatmaps — RUN IN CURSOR / RStudio.
# Replaces every gseapy.prerank result with genuine fgseaMultilevel
# (nPermSimple=10000, eps=0, scoreType="std", seed=42), on ONE identical
# Hallmark + KEGG + Reactome + GO:BP universe so the whole unbiased set is
# comparable AND consistent with the focused 58-panel (which also uses GO:BP).
#
# Produces three tidy CSVs (columns: contrast_id, pathway, collection, ES, NES,
# pval, padj[BH per contrast], size, leadingEdge [+ category for the focused set]):
#
#   1. data/MSigDB_HKR_v3_MULTILEVEL_topupdown_7contrasts.csv
#        universe: H+KEGG+Reactome+GO:BP (minSize 15, maxSize 500)
#        contrasts: E4-vs-E3; DAB1_KD & DAB2_KD (E3 / E4 / pooled)
#        -> single-column top-15-up/10-down heatmaps
#
#   2. data/MSigDB_HKR_v3_MULTILEVEL_3col_persiRNA.csv
#        universe: H+KEGG+Reactome+GO:BP (minSize 15, maxSize 500)
#        contrasts: per-siRNA E4-KD / E3-KD / interaction for 893, 96, 98
#        -> 3-column (APOE4 KD | APOE3 KD | interaction) heatmaps, now WITH KEGG
#
#   3. data/FOCUSED_v3_MULTILEVEL_7contrasts.csv
#        universe: the curated 58-pathway MG panel (H+KEGG+Reactome+GO:BP;
#        minSize 5, maxSize 2000), categories from FOCUSED_panel_definition_58.csv
#        -> the MG-keyword MAIN focused categorized heatmap
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
ARCH      <- "/Volumes/VERBATIM HD/FINAL FIGURES & SCRIPTS_JUNE2026/RELN_v2_redo_Wald/_archived_v3_strict80pct_filter"
RANK_ARCH <- file.path(ARCH, "ranked_Wald_v3_E3E4")     # pooled / genotype-split ranks
RANK_MAIN <- file.path(PROJ, "ranked_Wald")             # per-siRNA ranks (893/96/98)
OUT_DIR   <- file.path(PROJ, "Revised Figures_August2026/Unbiased_HKR_v3_MULTILEVEL/data")
PANEL_CSV <- file.path(OUT_DIR, "FOCUSED_panel_definition_58.csv")   # bundled panel def
# ---------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- gene sets ------------------------------------------------------------
cat("[1] loading MSigDB collections via msigdbr ...\n")
grab <- function(cat_, want=NULL){
  x <- msigdbr(species="Homo sapiens", category=cat_)
  sub <- if ("gs_subcat" %in% names(x)) x$gs_subcat else x$gs_subcollection
  if (!is.null(want)) x <- x[sub %in% want | grepl(want[1], sub), ]
  x
}
h    <- msigdbr(species="Homo sapiens", category="H")
c2   <- msigdbr(species="Homo sapiens", category="C2")
csub <- if ("gs_subcat" %in% names(c2)) c2$gs_subcat else c2$gs_subcollection
reac <- c2[csub == "CP:REACTOME", ]
kegg <- c2[grepl("^CP:KEGG", csub) & !grepl("MEDICUS", csub), ]           # KEGG legacy
gobp <- msigdbr(species="Homo sapiens", category="C5")
gsub <- if ("gs_subcat" %in% names(gobp)) gobp$gs_subcat else gobp$gs_subcollection
gobp <- gobp[gsub == "GO:BP", ]

mk <- function(df) split(toupper(df$gene_symbol), df$gs_name)
HKR   <- c(mk(h), mk(kegg), mk(reac))                    # unbiased universe
HKRGO <- c(HKR, mk(gobp))                                # + GO:BP for focused panel
cat(sprintf("    H=%d KEGG=%d Reactome=%d  -> HKR=%d sets  (+GO:BP=%d)\n",
            length(unique(h$gs_name)), length(unique(kegg$gs_name)),
            length(unique(reac$gs_name)), length(HKR), length(unique(gobp$gs_name))))

# ---- helpers --------------------------------------------------------------
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

# ---- SET 1: unbiased top-up/down (pooled + genotype), H+K+R+GO:BP --------
cat("[2] SET 1 — unbiased top-up/down (H+K+R+GO:BP) ...\n")
A <- c(BIONMG_E4_vs_E3_Ctrl="BIONMG_E4_vs_E3_Ctrl.rnk",
       BIONMG_DAB1_KD_vs_Ctrl_E3="BIONMG_DAB1_KD_vs_Ctrl_E3.rnk",
       BIONMG_DAB1_KD_vs_Ctrl_E4="BIONMG_DAB1_KD_vs_Ctrl_E4.rnk",
       BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4="BIONMG_DAB1_KD_vs_Ctrl_pooled_E3E4.rnk",
       BIONMG_DAB2_KD_vs_Ctrl_E3="BIONMG_DAB2_KD_vs_Ctrl_E3.rnk",
       BIONMG_DAB2_KD_vs_Ctrl_E4="BIONMG_DAB2_KD_vs_Ctrl_E4.rnk",
       BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4="BIONMG_DAB2_KD_vs_Ctrl_pooled_E3E4.rnk")
fwrite(run_many(A, RANK_ARCH, HKRGO, 15, 500),
       file.path(OUT_DIR, "MSigDB_HKR_v3_MULTILEVEL_topupdown_7contrasts.csv"))

# ---- SET 2: 3-column per-siRNA (E4-KD | E3-KD | interaction), H+K+R+GO:BP -
cat("[3] SET 2 — 3-column per-siRNA interaction (H+K+R+GO:BP) ...\n")
B <- c(E4_893_v3="BIONMG_893_vs_Ctrl_E4.rnk", E3_893_v3="BIONMG_893_vs_Ctrl_E3.rnk",
       int_893_v3="BIONMG_DAB1_response_E4_vs_E3_v3.rnk",   # coordinated DAB1 (893+894) interaction
       E4_894_v3="BIONMG_894_vs_Ctrl_E4.rnk", E3_894_v3="BIONMG_894_vs_Ctrl_E3.rnk",  # 894 shares int_893_v3
       E4_96_v3="BIONMG_96_vs_Ctrl_E4.rnk", E3_96_v3="BIONMG_96_vs_Ctrl_E3.rnk",
       int_96_v3="BIONMG_DAB2_response_E4_vs_E3_v3.rnk",
       E4_98_v3="BIONMG_98_vs_Ctrl_E4.rnk", E3_98_v3="BIONMG_98_vs_Ctrl_E3.rnk",
       int_98_v3="BIONMG_98_interaction_E4_vs_E3_v3.rnk")
fwrite(run_many(B, RANK_MAIN, HKRGO, 15, 500),
       file.path(OUT_DIR, "MSigDB_HKR_v3_MULTILEVEL_3col_persiRNA.csv"))

# ---- SET 3: curated 58-pathway MG focused panel (incl GO:BP) --------------
cat("[4] SET 3 — focused MG panel (H+K+R+GO:BP) ...\n")
panel <- fread(PANEL_CSV)                                  # pathway, category, setSize
focused_sets <- HKRGO[intersect(panel$pathway, names(HKRGO))]
missing <- setdiff(panel$pathway, names(HKRGO))
if (length(missing)) cat(sprintf("    NOTE: %d panel pathways not found in msigdbr: %s\n",
                                  length(missing), paste(missing, collapse=", ")))
focres <- run_many(A, RANK_ARCH, focused_sets, 5, 2000)
focres <- merge(focres, panel[, .(pathway, category)], by="pathway", all.x=TRUE)
fwrite(focres, file.path(OUT_DIR, "FOCUSED_v3_MULTILEVEL_7contrasts.csv"))

cat("\nDONE. Wrote 3 CSVs to", OUT_DIR, "\nThen run the three render_*.py scripts.\n")
