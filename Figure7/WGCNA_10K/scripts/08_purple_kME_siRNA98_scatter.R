#!/usr/bin/env Rscript
# =============================================================================
# WGCNA 10K — Purple module: kME × siRNA-98 gene-trait analysis
# =============================================================================
# For all 147 purple module genes, computes:
#   - kME: Pearson r with purple ME across all 39 samples (module hub score)
#   - Gene-trait r: Pearson r with siRNA-98 binary trait (1=siRNA-98, 0=Ctrl, n=16)
#   - BH-adjusted p-values
# Annotates RELN/APOE axis genes and AD risk genes.
# Produces scatter plot (kME vs gene-trait r) with labelled genes.
#
# Input:  /mnt/shared-workspace/datExpr_10K.rds
#         /mnt/shared-workspace/MEs_10K_p12.rds
#         /mnt/shared-workspace/moduleColors_10K_p12.rds
#         /mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv
#         /mnt/results/WGCNA/WGCNA_10K/tables/purple_module_DAB2_interactors.csv
#         <EDIT_INPUT_DIR>/MetaDataBatch.xlsx
#
# Output: /mnt/results/WGCNA/WGCNA_10K/tables/purple_kME_siRNA98_ranked_147genes.csv
#         /mnt/results/WGCNA/WGCNA_10K/figures/purple_kME_vs_siRNA98_scatter.png
# =============================================================================

library(readxl)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = FALSE)

# ---- Load data ----
datExpr      <- readRDS("/mnt/shared-workspace/datExpr_10K.rds")
MEs          <- readRDS("/mnt/shared-workspace/MEs_10K_p12.rds")
gene_assign  <- read.csv("/mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv")
purple_annot <- read.csv("/mnt/results/WGCNA/WGCNA_10K/tables/purple_module_DAB2_interactors.csv")
meta         <- read_xlsx("<EDIT_INPUT_DIR>/MetaDataBatch.xlsx")
meta$`Sample ID` <- as.integer(meta$`Sample ID`)

# Align metadata to MEs sample order
me_samples <- rownames(MEs)
meta_ord   <- meta[match(me_samples, as.character(meta$`Sample ID`)), ]
stopifnot(all(me_samples == as.character(meta_ord$`Sample ID`)))

# ---- Purple genes ----
purple_genes <- gene_assign[gene_assign$module == "purple", ]

# ---- kME ----
ME_purple  <- MEs[, "MEpurple"]
kme_purple <- sapply(purple_genes$ensembl_id, function(g)
  cor(datExpr[, g], ME_purple, use = "complete.obs"))

# ---- siRNA-98 gene-trait correlation ----
idx_ctrl  <- meta_ord$siRNA == "Ctrl"
idx_98    <- meta_ord$siRNA == "98"
valid_98  <- idx_ctrl | idx_98
trait_98  <- ifelse(idx_98, 1, ifelse(idx_ctrl, 0, NA))
n_98      <- sum(!is.na(trait_98))

corPvalueStudent <- function(r, n) {
  tstat <- r * sqrt(n - 2) / sqrt(1 - r^2)
  2 * pt(abs(tstat), df = n - 2, lower.tail = FALSE)
}

r_98    <- sapply(purple_genes$ensembl_id, function(g)
  cor(datExpr[valid_98, g], trait_98[valid_98], use = "complete.obs"))
p_98    <- sapply(r_98, corPvalueStudent, n = n_98)
padj_98 <- p.adjust(p_98, method = "BH")

# ---- Assemble results ----
axis_genes <- c("RELN","DAB2","DAB1","VLDLR","LRP8","LRP1","ITGB1","ITGB8","APOE")
results <- data.frame(
  ensembl_id     = purple_genes$ensembl_id,
  gene_symbol    = purple_genes$gene_symbol,
  kME_purple     = round(kme_purple, 4),
  r_siRNA98      = round(r_98, 4),
  pvalue_siRNA98 = signif(p_98, 4),
  padj_siRNA98   = signif(padj_98, 4),
  sig            = ifelse(padj_98 < 0.001, "***", ifelse(padj_98 < 0.01, "**",
                   ifelse(padj_98 < 0.05, "*", "ns"))),
  RELN_APOE_axis = purple_genes$gene_symbol %in% axis_genes
)
purple_sub <- purple_annot[, c("ensembl_id","AD_risk_gene","n_evidence",
  "STRING_DAB2_partner","DAB2_union","DAB1_union","LRP8_DAB1",
  "VLDLR_union","RELN_union","DAB2_lit_interactor")]
results <- merge(results, purple_sub, by = "ensembl_id", all.x = TRUE)
results <- results[order(results$r_siRNA98), ]

write.csv(results,
  "/mnt/results/WGCNA/WGCNA_10K/tables/purple_kME_siRNA98_ranked_147genes.csv",
  row.names = FALSE)

# ---- Scatter plot ----
plot_df <- results
plot_df$gene_symbol[is.na(plot_df$gene_symbol)] <- plot_df$ensembl_id[is.na(plot_df$gene_symbol)]

ad_risk_genes <- plot_df$gene_symbol[!is.na(plot_df$AD_risk_gene) & plot_df$AD_risk_gene]
top_r_genes   <- head(plot_df$gene_symbol[order(abs(plot_df$r_siRNA98), decreasing=TRUE)], 15)
label_genes   <- unique(c(plot_df$gene_symbol[plot_df$RELN_APOE_axis & !is.na(plot_df$RELN_APOE_axis)],
                           ad_risk_genes, top_r_genes))
plot_df$label <- ifelse(plot_df$gene_symbol %in% label_genes, plot_df$gene_symbol, NA)
plot_df$point_type <- factor(
  ifelse(plot_df$RELN_APOE_axis & !is.na(plot_df$RELN_APOE_axis), "RELN/APOE axis",
  ifelse(!is.na(plot_df$AD_risk_gene) & plot_df$AD_risk_gene, "AD risk gene",
  ifelse(plot_df$padj_siRNA98 < 0.05, "Significant (padj<0.05)", "Not significant"))),
  levels = c("RELN/APOE axis","AD risk gene","Significant (padj<0.05)","Not significant"))

colour_map <- c("RELN/APOE axis"="#E9ED4C","AD risk gene"="#FF9400",
                "Significant (padj<0.05)"="#9B59B6","Not significant"="grey70")
size_map   <- c("RELN/APOE axis"=4.5,"AD risk gene"=3.5,
                "Significant (padj<0.05)"=2.5,"Not significant"=2.0)

p <- ggplot(plot_df, aes(x=kME_purple, y=r_siRNA98, colour=point_type, size=point_type)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey60", linewidth=0.4) +
  geom_hline(yintercept=median(plot_df$r_siRNA98), linetype="dotted", colour="grey50", linewidth=0.3) +
  geom_vline(xintercept=median(plot_df$kME_purple), linetype="dotted", colour="grey50", linewidth=0.3) +
  geom_point(alpha=0.85) +
  geom_text_repel(aes(label=label), size=3.0,
    fontface=ifelse(plot_df$RELN_APOE_axis & !is.na(plot_df$RELN_APOE_axis), "bold.italic","plain"),
    colour="black", box.padding=0.4, point.padding=0.3, max.overlaps=40,
    segment.colour="grey50", segment.size=0.3, min.segment.length=0.2, na.rm=TRUE) +
  scale_colour_manual(values=colour_map, name=NULL) +
  scale_size_manual(values=size_map, name=NULL) +
  annotate("text", x=0.97, y=max(plot_df$r_siRNA98)*0.97,
    label="High hub\nUpregulated by DAB2 KD", hjust=1, size=2.8, colour="grey40", fontface="italic") +
  annotate("text", x=0.97, y=min(plot_df$r_siRNA98)*0.97,
    label="High hub\nDownregulated by DAB2 KD", hjust=1, size=2.8, colour="grey40", fontface="italic") +
  labs(title="Purple module: kME vs siRNA-98 gene-trait correlation",
    subtitle=paste0("147 purple genes  |  siRNA-98 vs Ctrl (n=16)  |  ",
      sum(plot_df$padj_siRNA98 < 0.05), "/147 significant (padj<0.05, BH)"),
    x="kME (correlation with purple ME, all 39 samples)",
    y="Gene-trait r  (siRNA-98 vs Ctrl)") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=12),
    plot.subtitle=element_text(size=9, colour="grey40"),
    legend.position="bottom", legend.text=element_text(size=9),
    panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="grey92")) +
  guides(colour=guide_legend(override.aes=list(size=3.5)), size="none")

ggsave("/mnt/results/WGCNA/WGCNA_10K/figures/purple_kME_vs_siRNA98_scatter.png",
  p, width=9, height=8, dpi=150, bg="white")
cat("Done.\n")

