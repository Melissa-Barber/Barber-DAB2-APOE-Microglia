# =============================================================================
# Volcano plots for DAB1 / DAB2 siRNA knockdown (vs control-siRNA), all
# genotypes pooled. Style matches the APOE4-vs-APOE3 volcano: blue = sig down,
# red = sig up (padj<0.05 & |log2FC|>1), grey = ns; dashed guides at |FC|=1
# and padj=0.05; knockdown target gene ringed and bold-labelled.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

deg_dir <- "/Volumes/VERBATIM HD/FINAL FIGURES & SCRIPTS_JUNE2026/DEG"
out_dir <- file.path(deg_dir, "DAB_KD_volcanoes"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LFC <- 1; AP <- 0.05; RED <- "#D6301F"; BLUE <- "#2C6FB2"; GREY <- "#B7B7B7"

volcano <- function(csv, target, title, subtitle = NULL) {
  d <- fread(csv)[!is.na(log2FoldChange) & !is.na(padj) & padj > 0]
  d[, nlp := -log10(padj)]
  d[, cls := fifelse(padj < AP & log2FoldChange >  LFC, "up",
              fifelse(padj < AP & log2FoldChange < -LFC, "down", "ns"))]
  xm <- max(3, ceiling(quantile(abs(d$log2FoldChange), 0.999)))
  lab <- rbind(head(d[cls == "up"][order(-nlp)], 8), head(d[cls == "down"][order(-nlp)], 8))
  lab <- lab[!is.na(hgnc_symbol) & hgnc_symbol != ""]
  tg  <- d[hgnc_symbol == target]

  ggplot(d, aes(log2FoldChange, nlp)) +
    geom_hline(yintercept = -log10(AP), linetype = "dashed", colour = "#9A9A9A", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "#9A9A9A", linewidth = 0.3) +
    geom_point(data = d[cls == "ns"],  colour = GREY, size = 0.7, alpha = 0.45) +
    geom_point(data = d[cls == "down"],colour = BLUE, size = 0.9, alpha = 0.8) +
    geom_point(data = d[cls == "up"],  colour = RED,  size = 0.9, alpha = 0.8) +
    geom_point(data = tg, shape = 21, colour = "black", fill = NA, size = 3, stroke = 0.9) +
    geom_text_repel(data = lab, aes(label = hgnc_symbol, colour = cls), fontface = "italic",
                    size = 2.6, seed = 3, max.overlaps = Inf, segment.size = 0.2, segment.colour = "#999") +
    geom_text_repel(data = tg, aes(label = hgnc_symbol), fontface = "bold.italic",
                    colour = "black", size = 3.0, seed = 3) +
    scale_colour_manual(values = c(up = RED, down = BLUE, ns = GREY), guide = "none") +
    coord_cartesian(xlim = c(-xm, xm)) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~fold~change~~(KD~vs~Ctrl)),
         y = expression(-log[10]~italic(p)[adj])) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11), plot.subtitle = element_text(size = 9))
}

specs <- list(
  list("Q3_DAB2_AllGenotypes_IsoformSpecific/Q3_DAB2_allGenotypes_96_vs_Ctrl.csv", "DAB2",
       "DAB2 siRNA-96 knockdown vs control", "targets DAB2-201 + DAB2-202", "volcano_DAB2_siRNA96.png"),
  list("Q3_DAB2_AllGenotypes_IsoformSpecific/Q3_DAB2_allGenotypes_98_vs_Ctrl.csv", "DAB2",
       "DAB2 siRNA-98 knockdown vs control", "targets DAB2-201", "volcano_DAB2_siRNA98.png"),
  list("Q3_DAB1_AllGenotypes_IsoformSpecific/Q3_DAB1_allGenotypes_893_vs_Ctrl.csv", "DAB1",
       "DAB1 siRNA-893 knockdown vs control", NULL, "volcano_DAB1_siRNA893.png"),
  list("Q3_DAB1_AllGenotypes_IsoformSpecific/Q3_DAB1_allGenotypes_894_vs_Ctrl.csv", "DAB1",
       "DAB1 siRNA-894 knockdown vs control", NULL, "volcano_DAB1_siRNA894.png"))

for (s in specs) {
  p <- volcano(file.path(deg_dir, s[[1]]), s[[2]], s[[3]], s[[4]])
  ggsave(file.path(out_dir, s[[5]]), p, width = 5.1, height = 3.5, dpi = 600, bg = "white")
  message("saved ", s[[5]])
}
