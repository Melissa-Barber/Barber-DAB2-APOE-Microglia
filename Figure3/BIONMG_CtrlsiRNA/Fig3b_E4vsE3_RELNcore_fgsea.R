## ============================================================================
## BIONMG_CtrlsiRNA_RELNcore_fgsea.R
## fgsea analysis: RELN-core-DAB1/2 gene sets in BIONMG APOE E4 vs E3 (Ctrl siRNA)
##
## Gene sets:
##   RELN_core_DAB1 (115 genes) = RELN_core_final_78genes ∪ RELN_DAB1_arm_optionA_37genes
##   RELN_core_DAB2 (113 genes) = RELN_core_final_78genes ∪ RELN_DAB2_arm_optionA_35genes
##
## Additional rows from v9 fgsea (pre-computed):
##   MG states (HM, DAM, HLA, CRM, IRM), Ribosomal (RM),
##   LDAM markers (Haney2024), APOE4 lipid iMG (Victor2022)
##
## Contrast: BIONMG E4 vs E3 (Ctrl siRNA)
##   Ranking: sign(log2FC) * -log10(pval)  [no Wald stat in this file]
##
## fgsea parameters: fgseaMultilevel, minSize=10, eps=0, nPermSimple=10000, set.seed(42)
##
## Outputs:
##   figures/BIONMG_APOE_E4vsE3_RELNcore_heatmap.png
##   data/BIONMG_E4vsE3_RELNcore_fgsea_results.csv
##
## Author: Melissa B. / Biomni session 2026-05-17
## ============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(ggplot2)
})

OUT_DIR      <- "."   # run from BIONMG_CtrlsiRNA/
UPLOAD_DIR   <- "<EDIT_INPUT_DIR>"
LISTS_DIR    <- file.path(UPLOAD_DIR, "RELN_lists")
V9_FGSEA     <- "<EDIT_PROJECT_ROOT>/fgsea_results/fgsea_all_contrasts_v9_clean.csv"
DESEQ2_FILE  <- "<EDIT_PROJECT_ROOT>/01_BIONMG_genotype_contrasts/BIONMG_E4vsE3_Ctrl_full_ranked.csv"

# ── 1. Gene sets ──────────────────────────────────────────────────────────
core <- read.csv(file.path(UPLOAD_DIR, "RELN_core_final_78genes.csv"))$gene_symbol
dab1 <- read.csv(file.path(LISTS_DIR,  "RELN_DAB1_arm_optionA_37genes.csv"))$gene_symbol
dab2 <- read.csv(file.path(LISTS_DIR,  "RELN_DAB2_arm_optionA_35genes.csv"))$gene_symbol
gene_sets <- list(
  RELN_core_DAB1 = union(core, dab1),
  RELN_core_DAB2 = union(core, dab2)
)
message("RELN_core_DAB1: ", length(gene_sets$RELN_core_DAB1), " genes")
message("RELN_core_DAB2: ", length(gene_sets$RELN_core_DAB2), " genes")

# ── 2. Ranked list ────────────────────────────────────────────────────────
df_b <- read.csv(DESEQ2_FILE, stringsAsFactors = FALSE)
ranked <- setNames(sign(df_b$log2FoldChange) * -log10(df_b$pvalue + 1e-300), df_b$symbol)
ranked <- ranked[!is.na(ranked) & !is.na(names(ranked)) & names(ranked) != ""]
ranked <- ranked[!duplicated(names(ranked))]
ranked <- sort(ranked, decreasing = TRUE)

# ── 3. fgsea ─────────────────────────────────────────────────────────────
set.seed(42)
res_reln <- fgseaMultilevel(gene_sets, ranked, minSize = 10, eps = 0, nPermSimple = 10000)
res_reln$contrast_id  <- "BIONMG_E4_vs_E3_Ctrl"
res_reln$leadingEdge  <- sapply(res_reln$leadingEdge, paste, collapse = ",")

# ── 4. v9 rows ────────────────────────────────────────────────────────────
v9 <- read.csv(V9_FGSEA, stringsAsFactors = FALSE)
v9_sets <- c("MG_HM","MG_DAM","MG_HLA","MG_CRM","MG_IRM",
             "Ribosomal response (RM)","Haney2024_LDAM_markers","Victor2022_APOE4_lipid_iMG")
v9_sub <- v9 %>%
  filter(pathway %in% v9_sets, contrast_id == "BIONMG_E4_vs_E3_Ctrl") %>%
  select(pathway, NES, padj, contrast_id)

# ── 5. Build plot data ────────────────────────────────────────────────────
row_order <- c("RELN_core_DAB1","RELN_core_DAB2",
               "MG_HM","MG_DAM","MG_HLA","MG_CRM","MG_IRM",
               "Ribosomal response (RM)","Haney2024_LDAM_markers","Victor2022_APOE4_lipid_iMG")
row_labels <- c(
  "RELN_core_DAB1"             = "RELN-core-DAB1",
  "RELN_core_DAB2"             = "RELN-core-DAB2",
  "MG_HM"                      = "MG Homeostatic (HM)",
  "MG_DAM"                     = "MG DAM",
  "MG_HLA"                     = "MG HLA",
  "MG_CRM"                     = "MG CRM",
  "MG_IRM"                     = "MG IRM",
  "Ribosomal response (RM)"    = "Ribosomal (RM)",
  "Haney2024_LDAM_markers"     = "LDAM markers (Haney)",
  "Victor2022_APOE4_lipid_iMG" = "APOE4 lipid iMG (Victor)"
)

plot_data <- bind_rows(v9_sub, res_reln %>% select(pathway, NES, padj, contrast_id)) %>%
  filter(pathway %in% row_order) %>%
  mutate(
    pathway_label = factor(row_labels[pathway], levels = rev(row_labels)),
    col_label     = "E4 vs E3",
    col_group     = "BIONMG\n(Ctrl siRNA)",
    sig_label     = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                              padj < 0.05  ~ "*",   TRUE ~ ""),
    NES_capped    = pmax(pmin(NES, 3), -3)
  )

# ── 6. Heatmap ────────────────────────────────────────────────────────────
p <- ggplot(plot_data, aes(x = col_label, y = pathway_label, fill = NES_capped)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sig_label), size = 4.0, vjust = 0.75, fontface = "bold") +
  geom_hline(yintercept = 8.5, color = "grey40", linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = 3.5, color = "grey40", linewidth = 0.6, linetype = "dashed") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-3, 3), name = "NES",
    breaks = c(-3, -1.5, 0, 1.5, 3),
    labels = c("\u2264-3", "-1.5", "0", "1.5", "\u22653")
  ) +
  facet_grid(. ~ col_group, scales = "free_x", space = "free_x") +
  labs(title = "RELN-core gene sets \u2014 APOE E4 vs E3 (BIONMG)", x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    text              = element_text(family = "Liberation Sans"),
    plot.title        = element_text(size = 10.5, face = "bold", hjust = 0.5,
                                     margin = margin(b = 8)),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 9.5),
    axis.text.y       = element_text(size = 9.5),
    axis.line         = element_blank(),
    axis.ticks        = element_blank(),
    strip.background  = element_rect(fill = "grey92", color = "grey70"),
    strip.text        = element_text(size = 9, face = "bold"),
    legend.position   = "right",
    legend.key.height = unit(1.1, "cm"),
    legend.key.width  = unit(0.35, "cm"),
    legend.title      = element_text(size = 8.5),
    legend.text       = element_text(size = 8),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

# ── 7. Save outputs ───────────────────────────────────────────────────────
ggsave(file.path(OUT_DIR, "figures", "BIONMG_APOE_E4vsE3_RELNcore_heatmap.png"),
       p, width = 4.0, height = 5.5, dpi = 300, bg = "white")

fgsea_out <- bind_rows(
  res_reln %>% select(pathway, NES, padj, size, leadingEdge, contrast_id),
  v9_sub %>% mutate(size = NA_integer_, leadingEdge = NA_character_)
) %>% arrange(match(pathway, row_order))
write.csv(fgsea_out,
          file.path(OUT_DIR, "data", "BIONMG_E4vsE3_RELNcore_fgsea_results.csv"),
          row.names = FALSE)

message("Done.")

