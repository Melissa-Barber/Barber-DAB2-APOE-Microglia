
# =============================================================================
# WGCNA 10K — Step 3: blockwiseModules (power = 12)
# =============================================================================
# Runs blockwiseModules() on the 10K gene expression matrix.
#
# IMPORTANT: This script must be run WITHOUT Bioconductor packages loaded
# (org.Hs.eg.db, clusterProfiler, etc.). Loading those packages overrides
# WGCNA's cor() with a stats::cor S4 generic that lacks the weights/cosine
# arguments WGCNA uses internally, causing a fatal error in the kME step.
# Run this script via: Rscript 03_blockwiseModules.R
#
# Parameters:
#   power          = 12   (R² = 0.846; chosen to reduce grey gene fraction)
#   networkType    = "signed hybrid"
#   TOMType        = "signed"
#   minModuleSize  = 30
#   mergeCutHeight = 0.25
#   deepSplit      = 2
#   corType        = "pearson"
#   maxBlockSize   = 10000  (single block)
#
# Result: 14 modules, 36.8% grey (3,684/10,000 genes)
#
# Input:  /mnt/shared-workspace/datExpr_10K.rds
# Output: /mnt/shared-workspace/net_10K_p12.rds
#         /mnt/shared-workspace/moduleColors_10K_p12.rds
#         /mnt/shared-workspace/MEs_10K_p12.rds
# =============================================================================

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)

datExpr_10K <- readRDS("/mnt/shared-workspace/datExpr_10K.rds")
cat("Matrix loaded:", nrow(datExpr_10K), "samples x", ncol(datExpr_10K), "genes\n")

set.seed(42)
net <- blockwiseModules(
  datExpr_10K,
  power             = 12,
  networkType       = "signed hybrid",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  deepSplit         = 2,
  corType           = "pearson",
  maxBlockSize      = 10000,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

# ---- Module summary ----
module_table <- sort(table(net$colors), decreasing = TRUE)
cat("\n--- Module summary (power=12) ---\n")
print(module_table)
cat("\nTotal modules (excl. grey):", sum(names(module_table) != "grey"), "\n")
cat("Genes in grey:", module_table["grey"],
    sprintf("(%.1f%%)\n", 100 * module_table["grey"] / sum(module_table)))

# ---- Save ----
saveRDS(net,          "/mnt/shared-workspace/net_10K_p12.rds")
saveRDS(net$colors,   "/mnt/shared-workspace/moduleColors_10K_p12.rds")
saveRDS(net$MEs,      "/mnt/shared-workspace/MEs_10K_p12.rds")
cat("Saved: net_10K_p12.rds, moduleColors_10K_p12.rds, MEs_10K_p12.rds\n")
