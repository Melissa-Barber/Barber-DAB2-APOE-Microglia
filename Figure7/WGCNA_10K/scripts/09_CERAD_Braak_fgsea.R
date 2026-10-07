# =============================================================================
# Script 09: CERAD/Braak fgsea + Updated AD Atlas Heatmap
# WGCNA 10K pipeline (power=12)
# =============================================================================
# Purpose:
#   Run fgsea for CERAD high vs low and Braak high vs low across 4 brain
#   regions (EC, ITG, PFC, V1) using pre-computed DESeq2 results.
#   Combine with existing PathAdv vs Early fgsea results.
#   Produce updated heatmap showing all 3 contrasts × 4 regions.
#
# Inputs:
#   /mnt/shared-workspace/deseq2_results/DESeq2_MG_{EC,ITG,PFC,V1}_{CERAD,Braak}_high_vs_low.csv
#   /mnt/results/WGCNA/WGCNA_10K/tables/AD_atlas_fgsea_10K_p12.csv  (PathAdv results)
#   /mnt/shared-workspace/module_gene_assignments_10K_p12.csv
#
# Outputs:
#   tables/AD_atlas_fgsea_10K_p12_withCERAD_Braak.csv  (160 rows, 3 contrasts)
#   figures/AD_atlas_fgsea_heatmap_10K_p12_withCERAD_Braak.png
# =============================================================================

library(fgsea)
library(data.table)
library(ggplot2)
options(stringsAsFactors = FALSE)

# ---- Paths ----
deg_dir    <- "/mnt/shared-workspace/deseq2_results"
results_dir <- "/mnt/results/WGCNA/WGCNA_10K"
dir.create(file.path(results_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# ---- Load module gene assignments ----
gene_assign <- read.csv("/mnt/shared-workspace/module_gene_assignments_10K_p12.csv",
                        stringsAsFactors = FALSE)
# Build gene sets (non-grey modules only)
all_modules <- unique(gene_assign$module)
all_modules <- all_modules[all_modules != "grey"]

module_genesets <- lapply(all_modules, function(m) {
  gene_assign$gene_symbol[gene_assign$module == m]
})
names(module_genesets) <- all_modules

cat("Gene sets loaded:", length(module_genesets), "modules\n")
cat("Module sizes:", paste(sapply(module_genesets, length), collapse=", "), "\n")

# ---- fgsea helper function ----
run_fgsea_contrast <- function(deg_file, contrast_name, region) {
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)
  
  # Rank metric: sign(log2FC) * -log10(pvalue)
  deg <- deg[!is.na(deg$pvalue) & deg$baseMean > 0, ]
  deg$rank_metric <- sign(deg$log2FoldChange) * (-log10(deg$pvalue))
  
  # Remove duplicates (keep highest |rank|)
  deg <- deg[order(abs(deg$rank_metric), decreasing = TRUE), ]
  deg <- deg[!duplicated(deg$gene), ]
  
  # Named rank vector
  ranks <- setNames(deg$rank_metric, deg$gene)
  
  # Run fgsea
  set.seed(42)
  res <- fgsea(
    pathways   = module_genesets,
    stats      = ranks,
    nPermSimple = 10000,
    minSize    = 10,
    maxSize    = Inf
  )
  
  res$region   <- region
  res$contrast <- contrast_name
  return(as.data.frame(res))
}

# ---- Run fgsea for all 8 CERAD/Braak contrasts ----
regions   <- c("EC", "ITG", "PFC", "V1")
contrasts <- c("CERAD_high_vs_low", "Braak_high_vs_low")

new_results <- list()
for (ct in contrasts) {
  for (region in regions) {
    key      <- paste0(region, "_", ct)
    deg_file <- file.path(deg_dir, paste0("DESeq2_MG_", region, "_", ct, ".csv"))
    
    if (!file.exists(deg_file)) {
      cat("WARNING: File not found:", deg_file, "\n")
      next
    }
    
    cat("Running fgsea:", key, "\n")
    res <- run_fgsea_contrast(deg_file, ct, region)
    new_results[[key]] <- res
  }
}

new_df <- do.call(rbind, new_results)
cat("New fgsea results:", nrow(new_df), "rows\n")

# ---- BH correction within each contrast ----
new_df$padj <- ave(new_df$pval, new_df$contrast, FUN = function(x) p.adjust(x, method = "BH"))

# ---- Load existing PathAdv results ----
pathadv <- read.csv(file.path(results_dir, "tables/AD_atlas_fgsea_10K_p12.csv"),
                    stringsAsFactors = FALSE)
cat("PathAdv results:", nrow(pathadv), "rows\n")

# ---- Align columns and combine ----
shared_cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size",
                 "leadingEdge", "region", "contrast")

# leadingEdge may be list in new_df — convert to character
if (is.list(new_df$leadingEdge)) {
  new_df$leadingEdge <- sapply(new_df$leadingEdge, paste, collapse=";")
}
if (is.list(pathadv$leadingEdge)) {
  pathadv$leadingEdge <- sapply(pathadv$leadingEdge, paste, collapse=";")
}

new_df_aligned    <- new_df[,    shared_cols]
pathadv_aligned   <- pathadv[,   shared_cols]
combined          <- rbind(pathadv_aligned, new_df_aligned)

# ---- Global BH correction across all 160 rows ----
combined$padj_global <- p.adjust(combined$pval, method = "BH")

cat("Combined:", nrow(combined), "rows\n")
cat("Significant (padj_global<0.05):", sum(combined$padj_global < 0.05, na.rm=TRUE), "\n")

# ---- Save combined CSV ----
write.csv(combined,
          file.path(results_dir, "tables/AD_atlas_fgsea_10K_p12_withCERAD_Braak.csv"),
          row.names = FALSE)
cat("Combined CSV saved.\n")

# ---- Generate heatmap ----
# Module order: by mean NES in PathAdv (most enriched → most depleted)
pathadv_nes <- tapply(combined$NES[combined$contrast == "PathAdv_vs_Early"],
                      combined$pathway[combined$contrast == "PathAdv_vs_Early"], mean, na.rm=TRUE)
mod_order   <- names(sort(pathadv_nes, decreasing = TRUE))
missing_mods <- setdiff(unique(combined$pathway), mod_order)
mod_order   <- c(mod_order, missing_mods)

contrast_order  <- c("PathAdv_vs_Early", "CERAD_high_vs_low", "Braak_high_vs_low")
contrast_labels <- c(
  "PathAdv_vs_Early"  = "Pathology\nAdv vs Early",
  "CERAD_high_vs_low" = "CERAD\nHigh vs Low",
  "Braak_high_vs_low" = "Braak\nHigh vs Low"
)
region_order <- c("EC", "ITG", "PFC", "V1")

full_grid <- expand.grid(
  pathway  = unique(combined$pathway),
  region   = region_order,
  contrast = contrast_order,
  stringsAsFactors = FALSE
)

plot_df <- merge(full_grid, combined[, c("pathway","region","contrast","NES","padj_global","size")],
                 by = c("pathway","region","contrast"), all.x = TRUE)

plot_df$sig_label <- ifelse(is.na(plot_df$padj_global), "",
                     ifelse(plot_df$padj_global < 0.001, "***",
                     ifelse(plot_df$padj_global < 0.01,  "**",
                     ifelse(plot_df$padj_global < 0.05,  "*", ""))))

# Flag turquoise CERAD/Braak cells (large gene set caveat)
plot_df$caveat <- plot_df$pathway == "turquoise" & plot_df$contrast != "PathAdv_vs_Early"

plot_df$pathway  <- factor(plot_df$pathway,  levels = rev(mod_order))
plot_df$region   <- factor(plot_df$region,   levels = region_order)
plot_df$contrast <- factor(plot_df$contrast, levels = contrast_order,
                           labels = unname(contrast_labels[contrast_order]))
plot_df$x_label  <- plot_df$region

mod_colours <- c(
  yellow="goldenrod2", red="firebrick", greenyellow="yellowgreen",
  green="forestgreen", purple="purple4", brown="saddlebrown",
  cyan="cyan4", black="black", magenta="magenta3", pink="hotpink",
  salmon="salmon", tan="tan3", blue="steelblue", turquoise="turquoise4"
)
axis_cols <- mod_colours[rev(mod_order)]
axis_cols[is.na(axis_cols)] <- "grey30"

plot_df$NES_capped <- pmax(pmin(plot_df$NES, 3), -3)

p <- ggplot(plot_df, aes(x = x_label, y = pathway, fill = NES_capped)) +
  geom_tile(aes(colour = caveat), linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 3.2, vjust = 0.75, na.rm = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-3, 3), na.value = "grey93",
    name = "NES"
  ) +
  scale_colour_manual(
    values = c("TRUE" = "grey50", "FALSE" = "white"),
    guide  = "none"
  ) +
  facet_grid(. ~ contrast, scales = "free_x", space = "free_x") +
  labs(
    title    = "WGCNA Module Enrichment in Human Microglia AD Progression",
    subtitle = paste0("fgsea (sign(log2FC)\u00d7\u2212log10p rank)  |  * padj_global<0.05  ** <0.01  *** <0.001\n",
                      "Grey border on turquoise (CERAD/Braak) = large gene set (n=1,468); interpret with caution"),
    x = "Brain region", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y      = element_text(size = 10, face = "bold", colour = axis_cols),
    axis.text.x      = element_text(size = 9),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "#ECE9E2", colour = NA),
    panel.grid       = element_blank(),
    panel.spacing    = unit(0.8, "lines"),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40", lineheight = 1.3),
    legend.position  = "right"
  )

ggsave(file.path(results_dir, "figures/AD_atlas_fgsea_heatmap_10K_p12_withCERAD_Braak.png"),
       p, width = 14, height = 7, dpi = 150, bg = "white")
cat("Heatmap saved.\n")

# ---- Print summary ----
cat("\n=== SUMMARY: Significant hits (padj_global < 0.05) ===\n")
sig <- combined[combined$padj_global < 0.05, ]
sig_summary <- sig[order(sig$padj_global), c("pathway","contrast","region","NES","padj_global")]
print(sig_summary, row.names = FALSE)
