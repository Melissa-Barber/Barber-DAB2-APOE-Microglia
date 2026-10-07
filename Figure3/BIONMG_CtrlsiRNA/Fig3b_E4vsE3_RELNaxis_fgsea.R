## ============================================================================
## BIONMG APOE e4 vs e3 — pathway-enrichment HEATMAP (WALD-ranked, KD-consistent)
##   input  : Wald-stat .rnk (gene <TAB> stat)         <- same as DAB1/DAB2 KD runs
##   engine : fgseaMultilevel, minSize=15, maxSize=500, eps=0, nPermSimple=10000, seed=42
##   sets   : Hallmark + KEGG + Reactome + GO:BP
##   output : all-black heatmap at 10.77 cm height + results table (with leading edges)
## Paste into R / Cursor and run.  Needs: msigdbr, fgsea, dplyr, readr, ggplot2
##   if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
##   BiocManager::install("fgsea"); install.packages(c("msigdbr","dplyr","readr","ggplot2"))
## ============================================================================

suppressPackageStartupMessages({
  library(msigdbr); library(fgsea); library(dplyr); library(readr); library(ggplot2)
})
set.seed(42)

## ---- 0. INPUT (edit if your drive path differs) ----------------------------
rnk_path <- "/Volumes/VERBATIM HD/Final Figures  June72026_v3heatmapsfinal/ranked_Wald/BIONMG_E4_vs_E3_Ctrl.rnk"
out_dir  <- "fgsea_E4vsE3_axis_heatmap_WALD"
dir.create(file.path(out_dir,"data"),    recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out_dir,"figures"), recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(rnk_path))

## ---- 1. Ranked vector: DESeq2 Wald stat (base R; dedup max|stat|, order) ----
d <- read.table(rnk_path, header=FALSE, sep="\t", quote="",
                stringsAsFactors=FALSE, col.names=c("gene","stat"))
d <- d[!is.na(d$gene) & !is.na(d$stat), ]
d$gene <- toupper(d$gene)
d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$gene), ]; d <- d[order(-d$stat), ]
ranked_vec <- setNames(d$stat, d$gene)

## ---- 2. Collections: Hallmark + KEGG + Reactome + GO:BP --------------------
h_df  <- msigdbr(species="Homo sapiens", collection="H")                                 %>% select(gs_name, gene_symbol)
kg_df <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:KEGG_LEGACY") %>% select(gs_name, gene_symbol)
re_df <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME")    %>% select(gs_name, gene_symbol)
bp_df <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP")           %>% select(gs_name, gene_symbol)
coll_df <- bind_rows(h_df, kg_df, re_df, bp_df) %>% mutate(gene_symbol = toupper(gene_symbol))
gs_list <- split(coll_df$gene_symbol, coll_df$gs_name)

## ---- 3. Multilevel fgsea ----------------------------------------------------
run_fgsea <- function(gs_df, label){
  res <- fgseaMultilevel(split(toupper(gs_df$gene_symbol), gs_df$gs_name), ranked_vec,
                         minSize=15, maxSize=500, eps=0, nPermSimple=10000)
  res$collection <- label; res
}
fgsea_all <- bind_rows(run_fgsea(h_df,"Hallmark"), run_fgsea(kg_df,"KEGG"),
                       run_fgsea(re_df,"Reactome"), run_fgsea(bp_df,"GO:BP"))
fgsea_all %>% mutate(leadingEdge = vapply(leadingEdge, paste, collapse="|", FUN.VALUE=character(1))) %>%
  write_csv(file.path(out_dir,"data","BIONMG_E4vsE3_fgsea_full_withLE_WALD.csv"))

## ---- 4. RELN/APOE-axis leading-edge genes (shown in each cell) --------------
AXIS_SHOW <- toupper(c("RELN","DAB1","DAB2","APOE","LRP8","VLDLR","LRP1","LRP2","ITGB1","ITGB8"))
fgsea_ann <- fgsea_all %>%
  mutate(shown    = lapply(leadingEdge, function(le) AXIS_SHOW[AXIS_SHOW %in% le]),
         k_axis   = lengths(shown),
         gene_lab = vapply(shown, function(g) if(length(g)) paste(g,collapse=", ") else "—", character(1)))

## ---- 5. Microglia keyword filter + exclude + Jaccard 0.5 collapse ----------
mg_kw <- paste0("immun|inflamm|cytokine|interferon|interleukin|nf.kb|nfkb|toll|tlr|innate|adaptive|",
  "complement|chemokine|mhc|antigen|phagocyt|lipid|cholesterol|lipoprotein|fatty.acid|",
  "triglyceride|apoe|lxr|sterol|microgli|myeloid|macrophage|monocyte|alzheimer|",
  "neurodegenerat|amyloid|tau|synapse|lysosom|autophagy|apoptosis|pi3k|mtor|mapk|",
  "jak.stat|nfat|trem|dap12|oxidative|reactive.oxygen|ros|stress|migration|adhesion|",
  "integrin|matrix.metalloprotein|ecm|extracellular.matrix|endocytosis|exocytosis|",
  "vesicle|phagosome|il[0-9]|tnf|ifn|ccl|cxcl|csf|tgf|neuroinflam|glia|astrocyte|",
  "oligodendro|neuron|brain|cns|proteasom|ubiquitin|protein.degradation|",
  "mitochondri|energy.metabol|glycolysis|oxidative.phosphorylation|cell.death|",
  "necrosis|pyroptosis|ferroptosis|dna.damage|genome.instability|senescence")
exclude <- c("GOBP_METANEPHROS_DEVELOPMENT","GOBP_HOMOPHILIC_CELL_CELL_ADHESION",
             "GOBP_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION")
mg <- fgsea_ann %>%
  filter(padj < 0.05, grepl(mg_kw, pathway, ignore.case=TRUE), !pathway %in% exclude) %>%
  arrange(padj)
jacc <- function(a,b) length(intersect(a,b))/length(union(a,b))
keep <- character(0); removed <- character(0)
for (p in mg$pathway){
  if (p %in% removed) next
  keep <- c(keep, p)
  for (q in mg$pathway){
    if (q==p || q %in% removed || q %in% keep) next
    if (jacc(gs_list[[p]], gs_list[[q]]) > 0.5) removed <- c(removed, q)
  }
}
mg_collapsed <- mg %>% filter(pathway %in% keep)

## ---- 6. Top 10 up + all significant down -----------------------------------
top_up   <- mg_collapsed %>% filter(NES>0) %>% arrange(padj) %>% head(10)
top_down <- mg_collapsed %>% filter(NES<0) %>% arrange(padj)

clean <- function(x){
  x <- gsub("^GOBP_|^REACTOME_|^HALLMARK_|^KEGG_(LEGACY_)?","",x)
  x <- tools::toTitleCase(tolower(gsub("_"," ", x)))
  x <- gsub("\\bEcm\\b","ECM", x);    x <- gsub("\\bTnfa\\b","TNFα", x)
  x <- gsub("\\bNfkb\\b","NF-κB", x); x <- gsub("\\bDna\\b","DNA", x)
  x <- gsub("\\bApoe\\b","APOE", x);  x <- gsub("\\bTlr\\b","TLR", x); x <- gsub("\\bMhc\\b","MHC", x)
  x
}
# order: up by NES ascending (top), then down by NES descending (below)
up   <- top_up   %>% arrange(NES)
down <- top_down %>% arrange(desc(NES))
panel <- bind_rows(up, down) %>%
  mutate(stars = case_when(padj<0.001~"***", padj<0.01~"**", padj<0.05~"*", TRUE~""),
         nes_lab  = sprintf("NES = %+.2f    %s", NES, stars),
         axis_lab = sprintf("%s   (%d/%d)", gene_lab, k_axis, size),
         name = clean(pathway))
N <- nrow(panel); panel$yy <- N:1               # first row at top
sep_y <- N - nrow(up) + 0.5                       # dashed line between up / down

write_csv(panel %>% select(pathway, collection, NES, padj, size, stars, k_axis, gene_lab),
          file.path(out_dir,"data","E4vsE3_axis_heatmap_table.csv"))

## ---- 7. All-black heatmap, lightened cells, ~10.77 cm ----------------------
INK <- "#1a1a1a"
g <- ggplot(panel, aes(x=0, y=yy, fill=NES)) +
  geom_tile(width=0.9, height=0.9, colour="white", linewidth=0.5) +
  geom_hline(yintercept=sep_y, linetype="dashed", colour="#9aa0a6", linewidth=0.35) +
  geom_text(aes(y=yy+0.20, label=nes_lab),  colour=INK, fontface="bold",        size=2.5) +
  geom_text(aes(y=yy-0.21, label=axis_lab), colour=INK, fontface="italic",      size=2.2) +
  scale_fill_gradient2(low="#7EA6CF", mid="#FFFFFF", high="#D27984",   # lightened RdBu
                       midpoint=0, limits=c(-3,3), breaks=seq(-3,3,1),
                       na.value="#EDEDED", name="NES") +
  scale_y_continuous(breaks=panel$yy, labels=panel$name, expand=expansion(add=0.6)) +
  scale_x_continuous(breaks=NULL, expand=c(0,0)) +
  labs(x=NULL, y=NULL,
       title="BIONMG iPSC-MG Ctrl-siRNA — APOE-ε4 vs ε3 pathway enrichment",
       subtitle=paste0("Hallmark · KEGG · Reactome · GO:BP | multilevel fgsea ",
                       "(min 15, max 500) · DESeq2 Wald-stat ranking\n",
                       "italic = RELN/APOE-axis genes in leading edge (k of N = set size) · ",
                       "*** <0.001  ** <0.01  * <0.05")) +
  theme_minimal(base_size=9) +
  theme(plot.title=element_text(face="bold", size=9.5),
        plot.subtitle=element_text(size=6, colour="#6b7078"),
        axis.text.y=element_text(size=7, colour="black"),
        axis.text.x=element_blank(), axis.ticks=element_blank(),
        panel.grid=element_blank(), legend.position="right",
        legend.title=element_text(size=7), legend.text=element_text(size=6),
        plot.margin=margin(4,6,3,4))
ggsave(file.path(out_dir,"figures","E4vsE3_axis_heatmap_WALD.png"), g,
       width=19/2.54, height=10.77/2.54, dpi=600, bg="white")
cat("Done. Table + figure in", out_dir, "\n")
