# =============================================================================
# DAB2-knockdown log2FC concordance across APOE genotypes (E3 vs E4)
# Per-isoform scatter of DESeq2 log2FoldChange, replacing the underpowered
# per-gene Genotype:siRNA interaction test with a transparent concordance plot.
#
# Input : Q3_DAB2_Separate/{E3,E4}_{96,98}_vs_Ctrl.csv  (DESeq2 results)
# Output: DAB2_isoform_genotype_concordance.png  (600 dpi)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

deg_dir <- "/Volumes/VERBATIM HD/FINAL FIGURES & SCRIPTS_JUNE2026/DEG/Q3_DAB2_Separate"
out_dir <- "/Volumes/VERBATIM HD/FINAL FIGURES & SCRIPTS_JUNE2026/DEG/DAB2_isoform_genotype_scatter"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_res <- function(f) {
  d <- fread(file.path(deg_dir, f))
  d <- d[!is.na(log2FoldChange), .(ensembl_id, log2FoldChange, padj)]
  d
}

# Deming / orthogonal regression (equal error-variance assumption)
deming_slope <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  sxx <- mean((x-mx)^2); syy <- mean((y-my)^2); sxy <- mean((x-mx)*(y-my))
  slope <- (syy - sxx + sqrt((syy - sxx)^2 + 4*sxy^2)) / (2*sxy)
  c(slope = slope, intercept = my - slope*mx)
}

read_sym <- function(f) {
  d <- fread(file.path(deg_dir, f))
  d[!is.na(log2FoldChange), .(ensembl_id, log2FoldChange, padj, hgnc_symbol)]
}

make_panel <- function(si, title, nlab = 6) {
  e3 <- read_sym(sprintf("E3_%s_vs_Ctrl.csv", si)); setnames(e3, c("ensembl_id","l3","p3","s3"))
  e4 <- read_sym(sprintf("E4_%s_vs_Ctrl.csv", si)); setnames(e4, c("ensembl_id","l4","p4","s4"))
  m  <- merge(e3, e4, by = "ensembl_id")
  m[, sym := fifelse(!is.na(s3), s3, s4)]
  m[is.na(p3), p3 := 1]; m[is.na(p4), p4 := 1]
  m[, sig3 := p3 < 0.05]; m[, sig4 := p4 < 0.05]
  m[, grp := fifelse(sig3 & sig4, "both",
              fifelse(sig3 | sig4, "one", "ns"))]

  both <- m[grp == "both"]
  x <- both$l3; y <- both$l4
  pr <- cor(x, y);  sr <- cor(x, y, method = "spearman")
  dm <- deming_slope(x, y)
  conc <- mean(sign(x) == sign(y)) * 100
  lim <- min(max(abs(c(m$l3, m$l4)))*0.62, 6.3); lim <- max(lim, 3.2)

  stat <- sprintf("%s:   n = %d genes significant in both genotypes    ·    Pearson r = %.2f    ·    Spearman rho = %.2f    ·    %.0f%% sign-concordant    ·    Deming slope = %.2f",
                  title, nrow(both), pr, sr, conc, dm["slope"])

  # most genotype-divergent genes (largest |E4-E3|), named & annotated only
  dv <- both[!is.na(sym)]
  dv[, sym1 := tstrsplit(sym, ";", fixed = TRUE, keep = 1L)]
  dv <- dv[!grepl("^LOC", sym1)]
  dv[, adiff := abs(l4 - l3)]
  dv <- dv[order(-adiff)][seq_len(min(nlab, .N))]

  ggplot() +
    geom_point(data = m[grp=="ns"],   aes(l3, l4), colour = "#EDEDED", size = 0.5, alpha = 0.55) +
    geom_point(data = m[grp=="one"],  aes(l3, l4), colour = "#8C8C8C", size = 1.1, alpha = 0.9) +
    geom_hline(yintercept = 0, colour = "#D8D8D8", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#D8D8D8", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, colour = "#B0342A", linetype = "dashed", linewidth = 0.5) +
    geom_point(data = both, aes(l3, l4), colour = "#1F3B6E", size = 1.5, alpha = 0.8) +
    geom_abline(slope = dm["slope"], intercept = dm["intercept"], colour = "#1F3B6E", linewidth = 0.8) +
    geom_point(data = dv, aes(l3, l4), shape = 21, colour = "#B0342A", fill = NA, stroke = 0.7, size = 2.6) +
    geom_text_repel(data = dv, aes(l3, l4, label = sym1), fontface = "italic",
                    size = 2.7, colour = "#111111", segment.colour = "#888888",
                    segment.size = 0.3, min.segment.length = 0, max.overlaps = Inf,
                    box.padding = 0.5, seed = 7) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = title,
         x = expression(APOE3~~log[2]*FC~~(DAB2~KD~vs~Ctrl)),
         y = expression(APOE4~~log[2]*FC~~(DAB2~KD~vs~Ctrl))) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))

  list(plot = p, stat = stat)
}

rA <- make_panel("96", "siRNA-96 (DAB2-201 + DAB2-202)", nlab = 5)
rB <- make_panel("98", "siRNA-98 (DAB2-201)", nlab = 6)

# stats printed as a two-line caption BELOW the shared legend, not inside panels
cap <- paste(rA$stat, rB$stat, sep = "\n")

fig <- (rA$plot | rB$plot) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "DAB2-knockdown response is directionally concordant across APOE genotypes",
    caption = cap,
    theme   = theme(plot.title   = element_text(hjust = 0.5, size = 13),
                    plot.caption = element_text(hjust = 0.5, size = 8, lineheight = 1.4),
                    legend.position = "bottom")) &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "DAB2_isoform_genotype_concordance.png"),
       fig, width = 12.4, height = 6.2, dpi = 600, bg = "white")
message("Saved: ", file.path(out_dir, "DAB2_isoform_genotype_concordance.png"))
