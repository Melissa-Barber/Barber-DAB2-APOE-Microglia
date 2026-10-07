# =============================================================================
# Fig1_RELN_core_DAB1_DAB2_fgsea_heatmaps.R
#
# Purpose: Gene set enrichment analysis of RELN signalling pathway gene sets
#          (RELN_core+DAB1 and RELN_core+DAB2) across:
#            - Galatro 2017: aged vs young human microglia (Fig1a)
#            - AD Atlas (Wachter et al. 2024): PathAdv vs Early, CERAD high vs
#              low, and Braak high vs low in EC, ITG, PFC, V1 (Fig1b)
#
# Outputs:
#   - DESeq2_MG_{region}_{CERAD/Braak}_high_vs_low.csv  (pseudobulk DE results)
#   - fgsea_RELN_core_DAB1_DAB2_all_contrasts.csv        (combined fgsea results)
#   - Fig1a_Galatro_RELN_core_DAB1_DAB2.png
#   - Fig1b_{EC/ITG/PFC/V1}_RELN_core_DAB1_DAB2.png
#
# Author: [Author]
# Date:   [Date]
# =============================================================================


# =============================================================================
# 0. Libraries
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)          # sparse matrix I/O
  library(DESeq2)          # pseudobulk differential expression
  library(fgsea)           # gene set enrichment
  library(ComplexHeatmap)  # heatmap visualisation
  library(circlize)        # colorRamp2 for heatmap colour scale
  library(grid)            # unit() for heatmap dimensions
  library(dplyr)           # data wrangling
  library(readr)           # fast CSV I/O
})


# =============================================================================
# 1. Setup & Gene Set Construction
# =============================================================================

# ---- 1.1  Paths --------------------------------------------------------------

upload_dir   <- "<EDIT_INPUT_DIR>"
lists_dir    <- file.path(upload_dir, "RELN_lists")
atlas_dir    <- "/mnt/shared-workspace/ad_atlas_downloads"
out_dir      <- "/mnt/shared-workspace"
fig1a_dir    <- "/mnt/results/Figures/Fig1/Fig1a/figures"
fig1b_dir    <- "/mnt/results/Figures/Fig1/Fig1b/figures"

# Create output directories if they do not exist
dir.create(fig1a_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig1b_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1.2  Load gene lists ----------------------------------------------------

message("Loading RELN gene lists ...")

core_df <- read_csv(
  file.path(upload_dir, "RELN_core_final_78genes.csv"),
  show_col_types = FALSE
)
dab1_df <- read_csv(
  file.path(lists_dir, "RELN_DAB1_arm_optionA_37genes.csv"),
  show_col_types = FALSE
)
dab2_df <- read_csv(
  file.path(lists_dir, "RELN_DAB2_arm_optionA_35genes.csv"),
  show_col_types = FALSE
)

core_genes <- unique(core_df$gene_symbol)
dab1_genes <- unique(dab1_df$gene_symbol)
dab2_genes <- unique(dab2_df$gene_symbol)

message(sprintf("  Core genes : %d", length(core_genes)))
message(sprintf("  DAB1 genes : %d", length(dab1_genes)))
message(sprintf("  DAB2 genes : %d", length(dab2_genes)))

# ---- 1.3  Build combined gene sets by union ----------------------------------

RELN_core_DAB1 <- union(core_genes, dab1_genes)  # expected 115
RELN_core_DAB2 <- union(core_genes, dab2_genes)  # expected 113

message(sprintf("  RELN_core_DAB1 (union): %d genes", length(RELN_core_DAB1)))
message(sprintf("  RELN_core_DAB2 (union): %d genes", length(RELN_core_DAB2)))

# Named list for fgsea
gene_sets <- list(
  RELN_core_DAB1 = RELN_core_DAB1,
  RELN_core_DAB2 = RELN_core_DAB2
)


# =============================================================================
# 2. Pseudobulk DESeq2 — AD Atlas CERAD and Braak contrasts
# =============================================================================

regions <- c("EC", "ITG", "PFC", "V1")

# Helper: rank-deduplicate a named numeric vector (keep highest |stat| per gene)
dedup_stats <- function(stats_vec) {
  df <- data.frame(gene = names(stats_vec), stat = stats_vec,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat), ]
  df <- df[order(abs(df$stat), decreasing = TRUE), ]
  df <- df[!duplicated(df$gene), ]
  setNames(df$stat, df$gene)
}

# Storage for DESeq2 results (used later in fgsea)
deseq2_results <- list()

for (region in regions) {

  message(sprintf("\n--- Pseudobulk DESeq2: %s ---", region))

  # ---- 2.1  Read sparse matrix -----------------------------------------------

  mtx_path  <- file.path(atlas_dir, region, "matrix.mtx")
  row_path  <- file.path(atlas_dir, region, "row_annotation.txt")
  cell_path <- file.path(atlas_dir, region, "cell_annotation.txt")

  counts_mat <- readMM(mtx_path)                          # genes x cells
  row_ann    <- read.table(row_path,  header = FALSE,
                           stringsAsFactors = FALSE)
  cell_ann   <- read.table(cell_path, header = FALSE,
                           stringsAsFactors = FALSE)

  rownames(counts_mat) <- row_ann[, 1]   # gene symbols / IDs
  colnames(counts_mat) <- cell_ann[, 1]  # cell barcodes

  # ---- 2.2  Read per-region metadata -----------------------------------------

  meta_file <- file.path(
    upload_dir,
    sprintf("2026-05-17_Microglia_%s_metadata.csv", region)
  )
  meta <- read_csv(meta_file, show_col_types = FALSE)
  # Expected columns: X (barcode), Donor.ID, Braak, CERAD, Sex, Pathology.Stage
  colnames(meta)[colnames(meta) == "X"] <- "barcode"

  # ---- 2.3  Align cells -------------------------------------------------------

  shared_cells <- intersect(colnames(counts_mat), meta$barcode)
  message(sprintf("  Shared cells: %d / %d matrix, %d metadata",
                  length(shared_cells), ncol(counts_mat), nrow(meta)))

  counts_mat <- counts_mat[, shared_cells, drop = FALSE]
  meta        <- meta[match(shared_cells, meta$barcode), ]

  # ---- 2.4  Aggregate pseudobulk counts per donor ----------------------------

  donors <- unique(meta$Donor.ID)
  pb_list <- lapply(donors, function(d) {
    cells_d <- meta$barcode[meta$Donor.ID == d]
    Matrix::rowSums(counts_mat[, cells_d, drop = FALSE])
  })
  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- donors
  pb_mat <- as.matrix(pb_mat)

  # Donor-level metadata (one row per donor)
  donor_meta <- meta[!duplicated(meta$Donor.ID), ]
  donor_meta  <- donor_meta[match(donors, donor_meta$Donor.ID), ]
  rownames(donor_meta) <- donor_meta$Donor.ID

  # ---- 2.5  CERAD contrast: high vs low --------------------------------------
  #   high = frequent / moderate neuritic plaques
  #   low  = sparse / none

  message("  Running CERAD contrast ...")

  donor_meta$CERAD_group <- dplyr::case_when(
    tolower(donor_meta$CERAD) %in% c("frequent", "moderate") ~ "high",
    tolower(donor_meta$CERAD) %in% c("sparse", "none")       ~ "low",
    TRUE ~ NA_character_
  )

  cerad_donors <- donor_meta$Donor.ID[!is.na(donor_meta$CERAD_group)]
  cerad_meta   <- donor_meta[cerad_donors, , drop = FALSE]
  cerad_meta$CERAD_group <- factor(cerad_meta$CERAD_group,
                                   levels = c("low", "high"))
  cerad_meta$Sex <- factor(cerad_meta$Sex)

  cerad_counts <- pb_mat[, cerad_donors, drop = FALSE]

  dds_cerad <- DESeqDataSetFromMatrix(
    countData = cerad_counts,
    colData   = cerad_meta,
    design    = ~ Sex + CERAD_group
  )
  dds_cerad <- DESeq(dds_cerad, quiet = TRUE)
  res_cerad <- results(dds_cerad,
                       contrast = c("CERAD_group", "high", "low"),
                       independentFiltering = TRUE)
  res_cerad_df <- as.data.frame(res_cerad)
  res_cerad_df$gene <- rownames(res_cerad_df)

  cerad_out <- file.path(
    out_dir,
    sprintf("DESeq2_MG_%s_CERAD_high_vs_low.csv", region)
  )
  write_csv(res_cerad_df, cerad_out)
  message(sprintf("  Saved: %s", basename(cerad_out)))

  deseq2_results[[paste0(region, "_CERAD")]] <- res_cerad_df

  # ---- 2.6  Braak contrast: high (V/VI) vs low (0/I/II) ----------------------
  #   Exclude Braak III (intermediate) and "unk"

  message("  Running Braak contrast ...")

  donor_meta$Braak_group <- dplyr::case_when(
    as.character(donor_meta$Braak) %in% c("V", "VI")         ~ "high",
    as.character(donor_meta$Braak) %in% c("0", "I", "II")    ~ "low",
    TRUE ~ NA_character_   # III, unk, or other -> excluded
  )

  braak_donors <- donor_meta$Donor.ID[!is.na(donor_meta$Braak_group)]
  braak_meta   <- donor_meta[braak_donors, , drop = FALSE]
  braak_meta$Braak_group <- factor(braak_meta$Braak_group,
                                   levels = c("low", "high"))
  braak_meta$Sex <- factor(braak_meta$Sex)

  braak_counts <- pb_mat[, braak_donors, drop = FALSE]

  dds_braak <- DESeqDataSetFromMatrix(
    countData = braak_counts,
    colData   = braak_meta,
    design    = ~ Sex + Braak_group
  )
  dds_braak <- DESeq(dds_braak, quiet = TRUE)
  res_braak <- results(dds_braak,
                       contrast = c("Braak_group", "high", "low"),
                       independentFiltering = TRUE)
  res_braak_df <- as.data.frame(res_braak)
  res_braak_df$gene <- rownames(res_braak_df)

  braak_out <- file.path(
    out_dir,
    sprintf("DESeq2_MG_%s_Braak_high_vs_low.csv", region)
  )
  write_csv(res_braak_df, braak_out)
  message(sprintf("  Saved: %s", basename(braak_out)))

  deseq2_results[[paste0(region, "_Braak")]] <- res_braak_df
}


# =============================================================================
# 3. fgseaMultilevel
# =============================================================================

message("\n--- Running fgseaMultilevel ---")

# fgsea parameters
FGSEA_PARAMS <- list(
  minSize      = 10,
  eps          = 0,
  nPermSimple  = 10000,
  seed         = 42
)

# Helper: build ranked stats vector from a data frame
make_ranks <- function(df, gene_col, stat_col) {
  df <- df[!is.na(df[[stat_col]]), ]
  df <- df[order(abs(df[[stat_col]]), decreasing = TRUE), ]
  df <- df[!duplicated(df[[gene_col]]), ]
  stats_vec <- setNames(df[[stat_col]], df[[gene_col]])
  stats_vec
}

# Helper: run fgsea and annotate with contrast label
run_fgsea <- function(ranks, contrast_label) {
  set.seed(FGSEA_PARAMS$seed)
  res <- fgseaMultilevel(
    pathways    = gene_sets,
    stats       = ranks,
    minSize     = FGSEA_PARAMS$minSize,
    eps         = FGSEA_PARAMS$eps,
    nPermSimple = FGSEA_PARAMS$nPermSimple
  )
  res$contrast <- contrast_label
  as.data.frame(res)
}

all_fgsea <- list()

# ---- 3.1  Fig1a: Galatro aged vs young --------------------------------------

message("  Galatro aged_vs_young ...")
galatro_df <- read_csv(
  file.path(upload_dir, "DESeq2_Galatro_aged_vs_young.csv"),
  show_col_types = FALSE
)
galatro_ranks <- make_ranks(galatro_df, "gene_symbol", "rank_stat")
all_fgsea[["Galatro_aged_vs_young"]] <- run_fgsea(
  galatro_ranks, "Galatro_aged_vs_young"
)

# ---- 3.2  Fig1b: PathAdv_vs_Early for each region ---------------------------

for (region in regions) {
  message(sprintf("  PathAdv_vs_Early: %s ...", region))
  path_file <- file.path(
    upload_dir,
    sprintf("DESeq2_MG_%s_PathologyAdv_vs_Early.csv", region)
  )
  path_df    <- read_csv(path_file, show_col_types = FALSE)
  path_ranks <- make_ranks(path_df, "gene", "stat")
  label      <- sprintf("%s_PathAdv_vs_Early", region)
  all_fgsea[[label]] <- run_fgsea(path_ranks, label)
}

# ---- 3.3  Fig1b: CERAD and Braak contrasts (from pseudobulk DESeq2 above) ---

for (region in regions) {

  # CERAD
  message(sprintf("  CERAD high_vs_low: %s ...", region))
  cerad_df    <- deseq2_results[[paste0(region, "_CERAD")]]
  cerad_ranks <- make_ranks(cerad_df, "gene", "stat")
  label_c     <- sprintf("%s_CERAD_high_vs_low", region)
  all_fgsea[[label_c]] <- run_fgsea(cerad_ranks, label_c)

  # Braak
  message(sprintf("  Braak high_vs_low: %s ...", region))
  braak_df    <- deseq2_results[[paste0(region, "_Braak")]]
  braak_ranks <- make_ranks(braak_df, "gene", "stat")
  label_b     <- sprintf("%s_Braak_high_vs_low", region)
  all_fgsea[[label_b]] <- run_fgsea(braak_ranks, label_b)
}

# ---- 3.4  Combine and save --------------------------------------------------

fgsea_combined <- dplyr::bind_rows(all_fgsea)
# Convert leadingEdge list column to a semicolon-separated string for CSV
fgsea_combined$leadingEdge <- sapply(
  fgsea_combined$leadingEdge,
  function(x) paste(x, collapse = ";")
)

fgsea_out <- file.path(out_dir, "fgsea_RELN_core_DAB1_DAB2_all_contrasts.csv")
write_csv(fgsea_combined, fgsea_out)
message(sprintf("Saved combined fgsea results: %s", fgsea_out))


# =============================================================================
# 4. Heatmap Generation (ComplexHeatmap)
# =============================================================================

message("\n--- Generating heatmaps ---")

# ---- 4.1  Shared heatmap settings -------------------------------------------

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#D6604D"))

# Significance stars helper
sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    padj < 0.1   ~ ".",
    TRUE         ~ ""
  )
}

# Cell function: NES (bold, size 7) above centre; stars below centre
make_cell_fun <- function(nes_mat, padj_mat) {
  function(j, i, x, y, width, height, fill) {
    nes_val  <- nes_mat[i, j]
    padj_val <- padj_mat[i, j]
    stars    <- sig_stars(padj_val)

    # NES label (bold, fontsize 7) -- slightly above centre
    grid.text(
      label  = sprintf("%.2f", nes_val),
      x      = x,
      y      = y + unit(0.12, "cm"),
      gp     = gpar(fontsize = 7, fontface = "bold")
    )
    # Stars -- slightly below centre
    if (nchar(stars) > 0) {
      grid.text(
        label  = stars,
        x      = x,
        y      = y - unit(0.12, "cm"),
        gp     = gpar(fontsize = 7)
      )
    }
  }
}

# Row order: RELN_core_DAB1 first, then RELN_core_DAB2
row_order <- c("RELN_core_DAB1", "RELN_core_DAB2")

# Helper: extract NES and padj matrices for given contrasts
extract_matrices <- function(fgsea_df, contrast_labels) {
  nes_mat  <- matrix(NA_real_, nrow = 2, ncol = length(contrast_labels),
                     dimnames = list(row_order, contrast_labels))
  padj_mat <- matrix(1,        nrow = 2, ncol = length(contrast_labels),
                     dimnames = list(row_order, contrast_labels))

  for (gs in row_order) {
    for (ct in contrast_labels) {
      row_idx <- fgsea_df$pathway == gs & fgsea_df$contrast == ct
      if (any(row_idx)) {
        nes_mat[gs, ct]  <- fgsea_df$NES[row_idx][1]
        padj_mat[gs, ct] <- fgsea_df$padj[row_idx][1]
      }
    }
  }
  list(nes = nes_mat, padj = padj_mat)
}

# Helper: draw and save a heatmap
save_heatmap <- function(nes_mat, padj_mat, col_split, col_split_labels,
                         filename_base) {

  # Replace NA NES with 0 for display (padj stays 1 -> no stars)
  nes_display <- nes_mat
  nes_display[is.na(nes_display)] <- 0

  n_cols    <- ncol(nes_mat)
  ht_width  <- unit(2.2 * n_cols, "cm")
  ht_height <- unit(5, "cm")

  ht <- Heatmap(
    matrix            = nes_display,
    name              = "NES",
    col               = col_fun,
    cell_fun          = make_cell_fun(nes_display, padj_mat),
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    row_names_side    = "left",
    show_column_names = TRUE,
    column_names_rot  = 45,
    column_names_gp   = gpar(fontsize = 8),
    row_names_gp      = gpar(fontsize = 8),
    column_split      = col_split,
    column_title_gp   = gpar(fill = "grey90", fontsize = 8, border = "grey70"),
    width             = ht_width,
    height            = ht_height,
    heatmap_legend_param = list(
      title         = "NES",
      at            = c(-2.5, 0, 2.5),
      labels        = c("-2.5", "0", "2.5"),
      legend_height = unit(3, "cm")
    )
  )

  # Save to /mnt/shared-workspace/ first
  ws_path <- file.path(out_dir, paste0(filename_base, ".png"))
  png(ws_path, width = 2.2 * n_cols + 5, height = 8, units = "cm",
      res = 300, bg = "white")
  draw(ht, padding = unit(c(0.5, 0.5, 2.5, 2), "cm"))  # 2.5 cm bottom for rotated multi-line column labels
  dev.off()
  message(sprintf("  Saved: %s", ws_path))

  ws_path
}

# ---- 4.2  Fig1a: Galatro aged vs young (1 column) ---------------------------

message("  Fig1a: Galatro ...")

galatro_contrasts <- "Galatro_aged_vs_young"
mats_1a <- extract_matrices(fgsea_combined, galatro_contrasts)

col_split_1a <- factor(rep("Aged vs Young", length(galatro_contrasts)),
                       levels = "Aged vs Young")

ws_1a <- save_heatmap(
  nes_mat          = mats_1a$nes,
  padj_mat         = mats_1a$padj,
  col_split        = col_split_1a,
  col_split_labels = "Aged vs Young",
  filename_base    = "Fig1a_Galatro_RELN_core_DAB1_DAB2"
)

# Copy to results directory (use system cp — file.copy() produces 0-byte files on S3-backed /mnt/results/)
system(sprintf('cp "%s" "%s"', ws_1a, file.path(fig1a_dir, basename(ws_1a))))
message(sprintf("  Copied to: %s", fig1a_dir))

# ---- 4.3  Fig1b: one heatmap per region (3 columns each) --------------------

for (region in regions) {
  message(sprintf("  Fig1b: %s ...", region))

  contrasts_1b <- c(
    sprintf("%s_PathAdv_vs_Early",   region),
    sprintf("%s_CERAD_high_vs_low",  region),
    sprintf("%s_Braak_high_vs_low",  region)
  )

  mats_1b <- extract_matrices(fgsea_combined, contrasts_1b)

  # Column split: one label per column
  col_split_1b <- factor(
    c("Pathology stage", "CERAD", "Braak"),
    levels = c("Pathology stage", "CERAD", "Braak")
  )

  ws_1b <- save_heatmap(
    nes_mat          = mats_1b$nes,
    padj_mat         = mats_1b$padj,
    col_split        = col_split_1b,
    col_split_labels = c("Pathology stage", "CERAD", "Braak"),
    filename_base    = sprintf("Fig1b_%s_RELN_core_DAB1_DAB2", region)
  )

  # Copy to results directory (use system cp — file.copy() produces 0-byte files on S3-backed /mnt/results/)
  system(sprintf('cp "%s" "%s"', ws_1b, file.path(fig1b_dir, basename(ws_1b))))
  message(sprintf("  Copied to: %s", fig1b_dir))
}

message("\nAll done.")
