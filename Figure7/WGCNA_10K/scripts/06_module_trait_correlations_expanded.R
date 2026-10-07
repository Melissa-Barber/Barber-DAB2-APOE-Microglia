#!/usr/bin/env Rscript
# 06_module_trait_correlations_expanded.R
# Expanded module-trait correlations: individual siRNAs + E4 Ctrl-only comparison
# WGCNA 10K, power=12

library(readxl)
library(ggplot2)
library(reshape2)

# ---- Load data ----
MEs  <- readRDS("MEs_10K_p12.rds")
meta <- read_xlsx("MetaDataBatch.xlsx")
meta$`Sample ID` <- as.integer(meta$`Sample ID`)

# Match and order metadata to ME samples
me_samples  <- rownames(MEs)
meta_matched <- meta[meta$`Sample ID` %in% as.integer(me_samples), ]
meta_ord     <- meta_matched[match(me_samples, as.character(meta_matched$`Sample ID`)), ]
stopifnot(all(me_samples == as.character(meta_ord$`Sample ID`)))

# ---- Build trait matrix ----
idx_dab1 <- meta_ord$siRNA %in% c("893", "894")
idx_dab2 <- meta_ord$siRNA %in% c("96", "98")
idx_ctrl <- meta_ord$siRNA == "Ctrl"

traitData <- data.frame(
  DAB1siRNA = ifelse(idx_dab1, 1, ifelse(idx_ctrl, 0, NA)),
  DAB2siRNA = ifelse(idx_dab2, 1, ifelse(idx_ctrl, 0, NA)),
  E4_all    = ifelse(meta_ord$Genotype == "E4", 1, ifelse(meta_ord$Genotype == "E3", 0, NA)),
  siRNA_893 = ifelse(meta_ord$siRNA == "893", 1, ifelse(idx_ctrl, 0, NA)),
  siRNA_894 = ifelse(meta_ord$siRNA == "894", 1, ifelse(idx_ctrl, 0, NA)),
  siRNA_96  = ifelse(meta_ord$siRNA == "96",  1, ifelse(idx_ctrl, 0, NA)),
  siRNA_98  = ifelse(meta_ord$siRNA == "98",  1, ifelse(idx_ctrl, 0, NA)),
  E4_Ctrl   = ifelse(idx_ctrl & meta_ord$Genotype == "E4", 1,
                     ifelse(idx_ctrl & meta_ord$Genotype == "E3", 0, NA)),
  row.names = me_samples
)

# ---- Pearson correlations ----
corPvalueStudent <- function(cor, nSamples) {
  tstat <- cor * sqrt(nSamples - 2) / sqrt(1 - cor^2)
  2 * pt(abs(tstat), df = nSamples - 2, lower.tail = FALSE)
}

modules <- colnames(MEs)
traits  <- colnames(traitData)
cor_mat  <- matrix(NA, nrow=length(modules), ncol=length(traits), dimnames=list(modules, traits))
pval_mat <- matrix(NA, nrow=length(modules), ncol=length(traits), dimnames=list(modules, traits))
n_mat    <- matrix(NA, nrow=length(modules), ncol=length(traits), dimnames=list(modules, traits))

for (tr in traits) {
  trait_vec <- traitData[[tr]]
  valid     <- !is.na(trait_vec)
  n_valid   <- sum(valid)
  for (mod in modules) {
    r <- cor(MEs[[mod]][valid], trait_vec[valid], use="complete.obs")
    cor_mat[mod, tr]  <- r
    pval_mat[mod, tr] <- corPvalueStudent(r, n_valid)
    n_mat[mod, tr]    <- n_valid
  }
}
padj_mat <- apply(pval_mat, 2, function(p) p.adjust(p, method="BH"))
rownames(padj_mat) <- rownames(pval_mat)

# ---- Save table ----
results_long <- data.frame(
  module = rep(sub("^ME","",rownames(cor_mat)), ncol(cor_mat)),
  trait  = rep(colnames(cor_mat), each=nrow(cor_mat)),
  r      = as.vector(cor_mat),
  pvalue = as.vector(pval_mat),
  padj   = as.vector(padj_mat),
  n      = as.vector(n_mat)
)
results_long$sig <- ifelse(results_long$padj < 0.001, "***",
                    ifelse(results_long$padj < 0.01,  "**",
                    ifelse(results_long$padj < 0.05,  "*", "")))
write.csv(results_long, "module_trait_correlations_expanded_10K_p12.csv", row.names=FALSE)

# ---- Heatmap ----
clean_mod <- sub("^ME", "", rownames(cor_mat))
cor_df    <- as.data.frame(cor_mat); cor_df$module <- clean_mod
cor_long  <- melt(cor_df, id.vars="module", variable.name="trait", value.name="r")
padj_df   <- as.data.frame(padj_mat); padj_df$module <- clean_mod
padj_long <- melt(padj_df, id.vars="module", variable.name="trait", value.name="padj")
df <- merge(cor_long, padj_long, by=c("module","trait"))
df$sig_label <- ifelse(df$padj < 0.001, "***", ifelse(df$padj < 0.01, "**",
                ifelse(df$padj < 0.05, "*", "")))

mod_order <- c("purple","salmon","brown","green","greenyellow","pink",
               "blue","magenta","cyan","black","red","yellow","tan","turquoise","grey")
df$module <- factor(df$module, levels=rev(mod_order))

trait_labels <- c(
  DAB1siRNA="DAB1 siRNA\n(893+894 vs Ctrl)", DAB2siRNA="DAB2 siRNA\n(96+98 vs Ctrl)",
  E4_all="E4 vs E3\n(all samples)", siRNA_893="siRNA-893\nvs Ctrl",
  siRNA_894="siRNA-894\nvs Ctrl", siRNA_96="siRNA-96\nvs Ctrl",
  siRNA_98="siRNA-98\nvs Ctrl", E4_Ctrl="E4 vs E3\n(Ctrl only)"
)
df$trait <- factor(df$trait, levels=names(trait_labels), labels=unname(trait_labels))
df$facet <- ifelse(df$trait %in% unname(trait_labels)[1:3], "Combined traits",
            ifelse(df$trait %in% unname(trait_labels)[4:7], "Individual siRNAs", "Genotype (Ctrl)"))
df$facet <- factor(df$facet, levels=c("Combined traits","Individual siRNAs","Genotype (Ctrl)"))

mod_colours <- c(purple="purple", salmon="salmon", brown="saddlebrown", green="forestgreen",
  greenyellow="yellowgreen", pink="hotpink", blue="steelblue", magenta="magenta3",
  cyan="cyan4", black="black", red="firebrick", yellow="goldenrod2",
  tan="tan3", turquoise="turquoise4", grey="grey50")
axis_cols <- mod_colours[rev(mod_order)]

p <- ggplot(df, aes(x=trait, y=module, fill=r)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label=sig_label), size=3.5, vjust=0.75) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#D6604D",
    midpoint=0, limits=c(-1,1), name="Pearson r") +
  facet_grid(. ~ facet, scales="free_x", space="free_x") +
  labs(title="Module-Trait Correlations (WGCNA 10K, power=12)",
       subtitle="* padj<0.05  ** padj<0.01  *** padj<0.001  (BH correction per trait column)",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, size=8.5),
        axis.text.y=element_text(size=10, face="bold", colour=axis_cols),
        strip.text=element_text(size=10, face="bold"),
        strip.background=element_rect(fill="#ECE9E2", color=NA),
        panel.grid=element_blank(), legend.position="right",
        plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"),
        panel.spacing=unit(0.8,"lines"))

ggsave("module_trait_heatmap_expanded_10K_p12.png", p, width=14, height=7, dpi=150, bg="white")
cat("Done.\n")

