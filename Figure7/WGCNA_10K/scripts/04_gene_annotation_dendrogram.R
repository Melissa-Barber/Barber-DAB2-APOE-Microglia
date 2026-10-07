
# =============================================================================
# WGCNA 10K — Step 4: ENSG → Gene Symbol Annotation + Dendrogram
# =============================================================================
# Converts Ensembl IDs to HGNC gene symbols using org.Hs.eg.db, builds the
# module gene assignment table, and generates the gene dendrogram figure.
#
# Input:  /mnt/shared-workspace/moduleColors_10K_p12.rds
#         /mnt/shared-workspace/net_10K_p12.rds
# Output: /mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv
#         /mnt/results/WGCNA/WGCNA_10K/figures/dendrogram_10K_p12.png
# =============================================================================

library(WGCNA)
library(org.Hs.eg.db)
library(AnnotationDbi)
options(stringsAsFactors = FALSE)

moduleColors <- readRDS("/mnt/shared-workspace/moduleColors_10K_p12.rds")
net          <- readRDS("/mnt/shared-workspace/net_10K_p12.rds")

# ---- ENSG → Gene Symbol ----
ensg_ids   <- names(moduleColors)
ensg_clean <- sub("\\..*$", "", ensg_ids)   # strip version suffix (e.g. .1)

sym_map <- mapIds(
  org.Hs.eg.db,
  keys      = ensg_clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

cat("Mapped:", sum(!is.na(sym_map)), "/", length(sym_map), "genes to symbols\n")
cat("Unmapped:", sum(is.na(sym_map)), "\n")

# ---- Module assignment table ----
module_df <- data.frame(
  ensembl_id  = ensg_ids,
  gene_symbol = as.character(sym_map[ensg_clean]),
  module      = moduleColors,
  row.names   = NULL,
  stringsAsFactors = FALSE
)

dir.create("/mnt/results/WGCNA/WGCNA_10K/tables",  recursive = TRUE, showWarnings = FALSE)
write.csv(module_df,
          "/mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv",
          row.names = FALSE)
cat("Saved: module_gene_assignments_10K_p12.csv\n")

# Module size summary
cat("\nModule sizes:\n")
print(sort(table(moduleColors), decreasing = TRUE))

# ---- Dendrogram ----
dir.create("/mnt/results/WGCNA/WGCNA_10K/figures", recursive = TRUE, showWarnings = FALSE)
png("/mnt/results/WGCNA/WGCNA_10K/figures/dendrogram_10K_p12.png",
    width = 2400, height = 1200, res = 150)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene dendrogram and module colors (10K genes, power=12)"
)

dev.off()
cat("Saved: dendrogram_10K_p12.png\n")
