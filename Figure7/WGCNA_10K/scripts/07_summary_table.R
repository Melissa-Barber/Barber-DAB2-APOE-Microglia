# =============================================================================
# WGCNA 10K — Step 7: Master Summary Table
# =============================================================================
# Assembles a single master table combining:
#   - Module membership for all 10,000 genes (Ensembl ID, symbol, module)
#   - kME (module membership correlation) for each gene to its own module
#   - Gene–trait correlations for all 8 traits (siRNA-98, siRNA-96, siRNA-893,
#     siRNA-894, DAB1siRNA, DAB2siRNA, E4_all, E4_Ctrl)
#   - RELN/APOE axis gene flag
#   - AD risk gene flag and RELN pathway evidence (from purple interactor table)
#
# Input:  /mnt/shared-workspace/datExpr_10K.rds
#         /mnt/shared-workspace/MEs_10K_p12.rds
#         /mnt/shared-workspace/moduleColors_10K_p12.rds
#         /mnt/shared-workspace/coldata_39.rds
#         /mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv
#         /mnt/results/WGCNA/WGCNA_10K/tables/purple_module_DAB2_interactors.csv
#         <EDIT_INPUT_DIR>/MetaDataBatch.xlsx
#
# Output: /mnt/results/WGCNA/WGCNA_10K/tables/WGCNA_10K_master_gene_table_p12.csv
# =============================================================================

library(readxl)
options(stringsAsFactors = FALSE)

# ---- Load data ----
datExpr     <- readRDS("/mnt/shared-workspace/datExpr_10K.rds")
MEs         <- readRDS("/mnt/shared-workspace/MEs_10K_p12.rds")
modColors   <- readRDS("/mnt/shared-workspace/moduleColors_10K_p12.rds")
coldata     <- readRDS("/mnt/shared-workspace/coldata_39.rds")
gene_assign <- read.csv("/mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv",
                        stringsAsFactors = FALSE)
purple_annot <- read.csv("/mnt/results/WGCNA/WGCNA_10K/tables/purple_module_DAB2_interactors.csv",
                         stringsAsFactors = FALSE)
meta        <- read_xlsx("<EDIT_INPUT_DIR>/MetaDataBatch.xlsx")
meta$`Sample ID` <- as.integer(meta$`Sample ID`)

# ---- Align metadata to datExpr sample order ----
me_samples  <- rownames(MEs)
meta_ord    <- meta[match(me_samples, as.character(meta$`Sample ID`)), ]
stopifnot(all(me_samples == as.character(meta_ord$`Sample ID`)))

# ---- Build trait vectors (pairwise: siRNA vs Ctrl only) ----
idx_ctrl <- meta_ord$siRNA == "Ctrl"

trait_list <- list(
  DAB1siRNA = ifelse(meta_ord$siRNA %in% c("893","894"), 1, ifelse(idx_ctrl, 0, NA)),
  DAB2siRNA = ifelse(meta_ord$siRNA %in% c("96","98"),   1, ifelse(idx_ctrl, 0, NA)),
  E4_all    = ifelse(meta_ord$Genotype == "E4", 1, ifelse(meta_ord$Genotype == "E3", 0, NA)),
  siRNA_893 = ifelse(meta_ord$siRNA == "893", 1, ifelse(idx_ctrl, 0, NA)),
  siRNA_894 = ifelse(meta_ord$siRNA == "894", 1, ifelse(idx_ctrl, 0, NA)),
  siRNA_96  = ifelse(meta_ord$siRNA == "96",  1, ifelse(idx_ctrl, 0, NA)),
  siRNA_98  = ifelse(meta_ord$siRNA == "98",  1, ifelse(idx_ctrl, 0, NA)),
  E4_Ctrl   = ifelse(idx_ctrl & meta_ord$Genotype == "E4", 1,
                     ifelse(idx_ctrl & meta_ord$Genotype == "E3", 0, NA))
)

# ---- Compute kME: correlation of each gene with its own module eigengene ----
# Uses all 39 samples
cat("Computing kME for", ncol(datExpr), "genes...\n")
kme_vec <- setNames(numeric(ncol(datExpr)), colnames(datExpr))
for (gene in colnames(datExpr)) {
  mod <- modColors[gene]
  me_col <- paste0("ME", mod)
  if (me_col %in% colnames(MEs)) {
    kme_vec[gene] <- cor(datExpr[, gene], MEs[, me_col], use = "complete.obs")
  } else {
    kme_vec[gene] <- NA
  }
}
cat("kME computed.\n")

# ---- Compute gene-trait correlations for all 10K genes ----
corPvalueStudent <- function(r, n) {
  tstat <- r * sqrt(n - 2) / sqrt(1 - r^2)
  2 * pt(abs(tstat), df = n - 2, lower.tail = FALSE)
}

gene_trait_r    <- matrix(NA, nrow = ncol(datExpr), ncol = length(trait_list),
                           dimnames = list(colnames(datExpr), names(trait_list)))
gene_trait_pval <- matrix(NA, nrow = ncol(datExpr), ncol = length(trait_list),
                           dimnames = list(colnames(datExpr), names(trait_list)))

for (tr in names(trait_list)) {
  tv    <- trait_list[[tr]]
  valid <- !is.na(tv)
  n_v   <- sum(valid)
  for (gene in colnames(datExpr)) {
    r <- cor(datExpr[valid, gene], tv[valid], use = "complete.obs")
    gene_trait_r[gene, tr]    <- r
    gene_trait_pval[gene, tr] <- corPvalueStudent(r, n_v)
  }
}

# BH correction per trait column
gene_trait_padj <- apply(gene_trait_pval, 2, function(p) p.adjust(p, method = "BH"))
rownames(gene_trait_padj) <- rownames(gene_trait_pval)
cat("Gene-trait correlations computed.\n")

# ---- Assemble master table ----
master <- gene_assign  # ensembl_id, gene_symbol, module

# Add kME
master$kME_own_module <- kme_vec[master$ensembl_id]

# Add gene-trait r and padj
for (tr in names(trait_list)) {
  master[[paste0("r_", tr)]]    <- gene_trait_r[master$ensembl_id, tr]
  master[[paste0("padj_", tr)]] <- gene_trait_padj[master$ensembl_id, tr]
}

# Add RELN/APOE axis flag
axis_genes <- c("RELN","DAB2","DAB1","VLDLR","LRP8","LRP1","ITGB1","ITGB8","APOE")
master$RELN_APOE_axis <- master$gene_symbol %in% axis_genes

# Add AD risk and RELN pathway evidence from purple interactor table
purple_sub <- purple_annot[, c("ensembl_id","AD_risk_gene","n_evidence",
                                "STRING_DAB2_partner","DAB2_union","DAB1_union",
                                "LRP8_DAB1","VLDLR_union","RELN_union")]
master <- merge(master, purple_sub, by = "ensembl_id", all.x = TRUE)

# Sort by module then kME descending
master <- master[order(master$module, -abs(master$kME_own_module), na.last = TRUE), ]

# Save
write.csv(master,
          "/mnt/results/WGCNA/WGCNA_10K/tables/WGCNA_10K_master_gene_table_p12.csv",
          row.names = FALSE)
cat("Saved: WGCNA_10K_master_gene_table_p12.csv\n")
cat("Dimensions:", nrow(master), "x", ncol(master), "\n")

# ---- Print RELN/APOE axis gene summary ----
cat("\n=== RELN/APOE axis genes ===\n")
axis_rows <- master[master$RELN_APOE_axis == TRUE & !is.na(master$RELN_APOE_axis), ]
print(axis_rows[, c("gene_symbol","module","kME_own_module",
                    "r_DAB2siRNA","padj_DAB2siRNA","r_siRNA_98","padj_siRNA_98")],
      row.names = FALSE)
