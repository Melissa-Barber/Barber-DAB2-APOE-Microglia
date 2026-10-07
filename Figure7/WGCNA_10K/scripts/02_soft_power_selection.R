
# =============================================================================
# WGCNA 10K — Step 2: Soft Power Selection
# =============================================================================
# Runs pickSoftThreshold() on the 10K gene expression matrix to identify the
# optimal soft-thresholding power for a signed hybrid network.
# Target: scale-free topology R² ≥ 0.85.
#
# Selected power: 12 (R² = 0.846, mean connectivity = 6.4)
# Note: power=14 (R²=0.873) was also tested but produced 44.7% grey genes;
#       power=12 reduced grey to 36.8% with acceptable scale-free fit.
#
# Input:  /mnt/shared-workspace/datExpr_10K.rds
# Output: /mnt/shared-workspace/sft_10K.rds
#         /mnt/shared-workspace/selected_power_10K.rds
#         /mnt/results/WGCNA/WGCNA_10K/figures/soft_power_plot_10K.png
# =============================================================================

library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)

datExpr_10K <- readRDS("/mnt/shared-workspace/datExpr_10K.rds")
cat("Matrix loaded:", nrow(datExpr_10K), "samples x", ncol(datExpr_10K), "genes\n")

# ---- Soft power selection ----
powers <- 1:30
sft <- pickSoftThreshold(
  datExpr_10K,
  powerVector  = powers,
  networkType  = "signed hybrid",
  verbose      = 5
)

saveRDS(sft, "/mnt/shared-workspace/sft_10K.rds")

# Print results table
cat("\nSoft power results:\n")
print(sft$fitIndices[, c("Power", "SFT.R.sq", "mean.k.")])

# ---- Select power ----
# First power where R² >= 0.85; use 12 to balance grey fraction
selected_power <- 12
saveRDS(selected_power, "/mnt/shared-workspace/selected_power_10K.rds")
cat("\nSelected power:", selected_power, "\n")

# ---- Plot ----
fit_df <- sft$fitIndices
fit_df$signed_R2 <- -sign(fit_df$slope) * fit_df$SFT.R.sq

dir.create("/mnt/results/WGCNA/WGCNA_10K/figures", recursive = TRUE, showWarnings = FALSE)
png("/mnt/results/WGCNA/WGCNA_10K/figures/soft_power_plot_10K.png",
    width = 1600, height = 700, res = 150)
par(mfrow = c(1, 2))

# R² vs power
plot(fit_df$Power, fit_df$signed_R2, type = "n",
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit (R²)",
     main = "Scale-free topology fit")
text(fit_df$Power, fit_df$signed_R2, labels = fit_df$Power, cex = 0.9, col = "red")
abline(h = 0.85, col = "blue", lty = 2)
abline(v = selected_power, col = "darkgreen", lty = 2)

# Mean connectivity vs power
plot(fit_df$Power, fit_df$mean.k., type = "n",
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity")
text(fit_df$Power, fit_df$mean.k., labels = fit_df$Power, cex = 0.9, col = "red")
abline(v = selected_power, col = "darkgreen", lty = 2)

dev.off()
cat("Saved: soft_power_plot_10K.png\n")
