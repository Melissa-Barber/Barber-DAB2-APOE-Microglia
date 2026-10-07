
# =============================================================================
# WGCNA 10K — Step 1: Gene Selection by Variance
# =============================================================================
# Selects the top 10,000 most variable genes from the VST batch-corrected
# expression matrix (23,320 genes × 39 samples; E3+E4, no E2).
# Transposes to WGCNA format (samples × genes) and runs goodSamplesGenes QC.
#
# Input:  /mnt/shared-workspace/vst_batch_corrected_allGenes_noE2.rds
# Output: /mnt/shared-workspace/datExpr_10K.rds
#         /mnt/shared-workspace/top10k_genes.rds
#         /mnt/results/WGCNA/WGCNA_10K/tables/top10K_genes_by_variance.csv
# =============================================================================

library(WGCNA)
options(stringsAsFactors = FALSE)

# ---- Load VST matrix ----
vst <- readRDS("/mnt/shared-workspace/vst_batch_corrected_allGenes_noE2.rds")
cat("VST matrix dimensions:", nrow(vst), "genes x", ncol(vst), "samples\n")

# ---- Select top 10,000 genes by variance ----
gene_vars  <- apply(vst, 1, var)
top10k_idx <- order(gene_vars, decreasing = TRUE)[1:10000]
top10k_ids <- rownames(vst)[top10k_idx]

cat("Variance range of selected genes:",
    round(min(gene_vars[top10k_idx]), 3), "–",
    round(max(gene_vars[top10k_idx]), 3), "\n")

# ---- Transpose to WGCNA format: samples × genes ----
datExpr_10K <- t(vst[top10k_ids, ])
cat("datExpr_10K dimensions:", nrow(datExpr_10K), "samples x", ncol(datExpr_10K), "genes\n")

# ---- QC: check for bad samples/genes ----
gsg <- goodSamplesGenes(datExpr_10K, verbose = 3)
if (!gsg$allOK) {
  datExpr_10K <- datExpr_10K[gsg$goodSamples, gsg$goodGenes]
  cat("After QC:", nrow(datExpr_10K), "samples x", ncol(datExpr_10K), "genes\n")
} else {
  cat("All samples and genes passed QC.\n")
}

# ---- Save outputs ----
saveRDS(datExpr_10K, "/mnt/shared-workspace/datExpr_10K.rds")
saveRDS(top10k_ids,  "/mnt/shared-workspace/top10k_genes.rds")

dir.create("/mnt/results/WGCNA/WGCNA_10K/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(
  data.frame(ensembl_id = top10k_ids,
             variance   = gene_vars[top10k_ids],
             rank       = seq_along(top10k_ids)),
  "/mnt/results/WGCNA/WGCNA_10K/tables/top10K_genes_by_variance.csv",
  row.names = FALSE
)
cat("Saved datExpr_10K.rds, top10k_genes.rds, top10K_genes_by_variance.csv\n")
