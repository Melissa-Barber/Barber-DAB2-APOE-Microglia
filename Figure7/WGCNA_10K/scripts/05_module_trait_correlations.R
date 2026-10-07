
# =============================================================================
# WGCNA 10K — Step 5: Module–Trait Correlations
# =============================================================================
# Builds a binary trait matrix from sample metadata (DAB1_KD, DAB2_KD,
# APOE_E3, APOE_E4) and computes Pearson correlations between module
# eigengenes and traits. Generates a heatmap with significance stars.
#
# Traits:
#   DAB1_KD  — 1 if siRNA_class == "DAB1", else 0  (n=17)
#   DAB2_KD  — 1 if siRNA_class == "DAB2", else 0  (n=14)
#   APOE_E3  — 1 if Genotype == "E3", else 0        (n=22)
#   APOE_E4  — 1 if Genotype == "E4", else 0        (n=17)
#
# Note: APOE_E3 and APOE_E4 are perfectly collinear in this dataset
#       (E3 = 1 − E4), so they carry no independent information.
#
# Input:  /mnt/shared-workspace/MEs_10K_p12.rds
#         /mnt/shared-workspace/datExpr_10K.rds
#         <EDIT_INPUT_DIR>/DESeq2_input_coldata.csv
# Output: /mnt/results/WGCNA/WGCNA_10K/tables/module_trait_correlations_p12.csv
#         /mnt/results/WGCNA/WGCNA_10K/tables/module_trait_pvalues_p12.csv
#         /mnt/results/WGCNA/WGCNA_10K/figures/module_trait_heatmap_p12.png
#         /mnt/shared-workspace/moduleTraitCor_p12.rds
#         /mnt/shared-workspace/moduleTraitPval_p12.rds
#         /mnt/shared-workspace/traits_p12.rds
# =============================================================================

library(WGCNA)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = FALSE)

MEs         <- readRDS("/mnt/shared-workspace/MEs_10K_p12.rds")
datExpr_10K <- readRDS("/mnt/shared-workspace/datExpr_10K.rds")

# ---- Build trait matrix ----
coldata     <- read.csv("<EDIT_INPUT_DIR>/DESeq2_input_coldata.csv", row.names = 1)
coldata_sub <- coldata[rownames(coldata) %in% rownames(datExpr_10K), ]
coldata_sub <- coldata_sub[match(rownames(MEs), rownames(coldata_sub)), ]

traits <- data.frame(
  DAB1_KD = as.integer(coldata_sub$siRNA_class == "DAB1"),
  DAB2_KD = as.integer(coldata_sub$siRNA_class == "DAB2"),
  APOE_E3 = as.integer(coldata_sub$Genotype == "E3"),
  APOE_E4 = as.integer(coldata_sub$Genotype == "E4"),
  row.names = rownames(coldata_sub)
)
cat("Trait sums:", colSums(traits), "\n")

# ---- Correlations ----
MEs_noGrey      <- MEs[, !grepl("grey", colnames(MEs))]
nSamples        <- nrow(MEs_noGrey)
moduleTraitCor  <- cor(MEs_noGrey, traits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nSamples)

# Save tables
dir.create("/mnt/results/WGCNA/WGCNA_10K/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(moduleTraitCor),
          "/mnt/results/WGCNA/WGCNA_10K/tables/module_trait_correlations_p12.csv",
          row.names = TRUE)
write.csv(as.data.frame(moduleTraitPval),
          "/mnt/results/WGCNA/WGCNA_10K/tables/module_trait_pvalues_p12.csv",
          row.names = TRUE)

saveRDS(moduleTraitCor,  "/mnt/shared-workspace/moduleTraitCor_p12.rds")
saveRDS(moduleTraitPval, "/mnt/shared-workspace/moduleTraitPval_p12.rds")
saveRDS(traits,          "/mnt/shared-workspace/traits_p12.rds")

# ---- Heatmap ----
rownames(moduleTraitCor)  <- sub("^ME", "", rownames(moduleTraitCor))
rownames(moduleTraitPval) <- sub("^ME", "", rownames(moduleTraitPval))

sig_stars <- ifelse(moduleTraitPval < 0.001, "***",
             ifelse(moduleTraitPval < 0.01,  "**",
             ifelse(moduleTraitPval < 0.05,  "*", "")))

cor_melt  <- reshape2::melt(moduleTraitCor,  varnames = c("Module", "Trait"), value.name = "Correlation")
star_melt <- reshape2::melt(sig_stars,       varnames = c("Module", "Trait"), value.name = "Stars")
pval_melt <- reshape2::melt(moduleTraitPval, varnames = c("Module", "Trait"), value.name = "Pvalue")
plot_df   <- merge(merge(cor_melt, star_melt, by = c("Module", "Trait")),
                   pval_melt, by = c("Module", "Trait"))

mod_order      <- rownames(moduleTraitCor)[order(moduleTraitCor[, "DAB2_KD"])]
plot_df$Module <- factor(plot_df$Module, levels = mod_order)
plot_df$Trait  <- factor(plot_df$Trait,
                         levels = c("DAB1_KD", "DAB2_KD", "APOE_E3", "APOE_E4"))
plot_df$label  <- paste0(sprintf("%.2f", plot_df$Correlation),
                         ifelse(plot_df$Stars != "", paste0("\n", plot_df$Stars), ""))

p <- ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.2, lineheight = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  scale_x_discrete(labels = c("DAB1_KD" = "DAB1 KD", "DAB2_KD" = "DAB2 KD",
                               "APOE_E3" = "APOE E3", "APOE_E4" = "APOE E4")) +
  labs(title    = "Module-Trait Correlations (WGCNA 10K, power=12)",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = "Module") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x     = element_text(angle = 30, hjust = 1, face = "bold"),
        axis.text.y     = element_text(face = "bold"),
        panel.grid      = element_blank(),
        plot.title      = element_text(face = "bold", size = 14),
        plot.subtitle   = element_text(size = 10, color = "grey40"),
        legend.position = "right")

dir.create("/mnt/results/WGCNA/WGCNA_10K/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("/mnt/results/WGCNA/WGCNA_10K/figures/module_trait_heatmap_p12.png",
       p, width = 7, height = 8, dpi = 150)
cat("Saved: module_trait_heatmap_p12.png\n")

# Print significant associations
cat("\nSignificant module-trait associations (p<0.05):\n")
sig_idx <- which(moduleTraitPval < 0.05, arr.ind = TRUE)
for (i in seq_len(nrow(sig_idx))) {
  mod   <- rownames(moduleTraitPval)[sig_idx[i, 1]]
  trait <- colnames(moduleTraitPval)[sig_idx[i, 2]]
  r     <- round(moduleTraitCor[sig_idx[i, 1], sig_idx[i, 2]], 3)
  p     <- signif(moduleTraitPval[sig_idx[i, 1], sig_idx[i, 2]], 3)
  cat(sprintf("  %-15s %-10s  r = %+.3f  p = %.4f\n", mod, trait, r, p))
}
