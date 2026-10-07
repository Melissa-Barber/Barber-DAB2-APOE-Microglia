## ============================================================================
## BIONMG APOE E4 vs E3 — Full Pathway fgsea (Microglia-Focused)
## Collections: MSigDB Hallmarks, GO:BP, Reactome (via msigdbr)
## Output: Ranked horizontal bar chart (top 10 up + significant down)
## ============================================================================

library(msigdbr); library(fgsea); library(dplyr); library(readr); library(ggplot2)

deseq_path <- "<EDIT_PROJECT_ROOT>/01_BIONMG_genotype_contrasts/BIONMG_E4vsE3_Ctrl_full_ranked.csv"
out_dir    <- "/mnt/results/Figures/Fig1/BIONMG_CtrlsiRNA"
dir.create(file.path(out_dir, "figures"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out_dir, "data"),    recursive=TRUE, showWarnings=FALSE)

# 1. Ranked vector: sign(log2FC) * -log10(pvalue)
deseq <- read_csv(deseq_path, show_col_types=FALSE)
ranked_vec <- deseq %>%
  filter(!is.na(symbol), !is.na(log2FoldChange), !is.na(pvalue), pvalue > 0) %>%
  mutate(m = sign(log2FoldChange) * -log10(pvalue)) %>%
  arrange(desc(m)) %>% distinct(symbol, .keep_all=TRUE) %>%
  { setNames(.$m, .$symbol) }

# 2. Gene sets
h_df  <- msigdbr(species="Homo sapiens", collection="H")                              %>% select(gs_name, gene_symbol)
bp_df <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP")      %>% select(gs_name, gene_symbol)
re_df <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME") %>% select(gs_name, gene_symbol)

# 3. fgsea per collection
run_fgsea <- function(gs, label) {
  set.seed(42)
  res <- fgseaMultilevel(split(gs$gene_symbol, gs$gs_name), ranked_vec,
                         minSize=15, maxSize=500, eps=0, nPermSimple=10000)
  res$collection <- label; res
}
fgsea_all <- bind_rows(run_fgsea(h_df,"Hallmark"), run_fgsea(bp_df,"GO:BP"), run_fgsea(re_df,"Reactome")) %>% select(-leadingEdge)
write_csv(fgsea_all, file.path(out_dir, "data", "BIONMG_E4vsE3_fgsea_full_results.csv"))

# 4. Microglia keyword filter
mg_kw <- paste0("immun|inflamm|cytokine|interferon|interleukin|nf.kb|nfkb|toll|tlr|innate|adaptive|",
  "complement|chemokine|mhc|antigen|phagocyt|lipid|cholesterol|lipoprotein|fatty.acid|",
  "triglyceride|apoe|lxr|sterol|microgli|myeloid|macrophage|monocyte|alzheimer|",
  "neurodegenerat|amyloid|tau|synapse|lysosom|autophagy|apoptosis|pi3k|mtor|mapk|",
  "jak.stat|nfat|trem|dap12|oxidative|reactive.oxygen|ros|stress|migration|adhesion|",
  "integrin|matrix.metalloprotein|ecm|extracellular.matrix|endocytosis|exocytosis|",
  "vesicle|phagosome|il[0-9]|tnf|ifn|ccl|cxcl|csf|tgf|neuroinflam|glia|astrocyte|",
  "oligodendro|neuron|brain|cns|proteasom|ubiquitin|protein.degradation|",
  "mitochondri|energy.metabol|glycolysis|oxidative.phosphorylation|cell.death|",
  "necrosis|pyroptosis|ferroptosis|dna.damage|genome.instability|senescence")

exclude <- c("GOBP_METANEPHROS_DEVELOPMENT","GOBP_HOMOPHILIC_CELL_CELL_ADHESION",
             "GOBP_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION")

mg_filtered <- fgsea_all %>%
  filter(padj < 0.05, grepl(mg_kw, pathway, ignore.case=TRUE), !pathway %in% exclude)

# 5. Jaccard redundancy collapse (threshold = 0.5)
gs_all <- split(bind_rows(h_df,bp_df,re_df)$gene_symbol, bind_rows(h_df,bp_df,re_df)$gs_name)
gs_sub <- gs_all[names(gs_all) %in% mg_filtered$pathway]
jaccard <- function(a,b) length(intersect(a,b))/length(union(a,b))
mg_sorted <- mg_filtered %>% arrange(padj)
keep <- character(0); removed <- character(0)
for (p in mg_sorted$pathway) {
  if (p %in% removed) next
  keep <- c(keep, p)
  if (!is.null(gs_sub[[p]])) for (q in mg_sorted$pathway) {
    if (q==p || q %in% removed || q %in% keep) next
    if (!is.null(gs_sub[[q]]) && jaccard(gs_sub[[p]],gs_sub[[q]]) > 0.5) removed <- c(removed,q)
  }
}
mg_collapsed <- mg_filtered %>% filter(pathway %in% keep)
write_csv(mg_collapsed, file.path(out_dir, "data", "BIONMG_E4vsE3_fgsea_MG_filtered.csv"))

# 6. Select top 10 up + all significant down
top_up   <- mg_collapsed %>% filter(NES > 0) %>% arrange(padj) %>% head(10)
top_down <- mg_collapsed %>% filter(NES < 0) %>% arrange(padj)
plot_data <- bind_rows(top_up, top_down)

# 7. Clean names and plot
clean_name <- function(x) {
  x <- gsub("^GOBP_|^REACTOME_|^HALLMARK_","",x); x <- gsub("_"," ",x)
  x <- tools::toTitleCase(tolower(x))
  for (w in c("Dna","Rna","Atp","Ldl","Hdl","Vldl","Apoe","Tnf","Tlr","Ecm"))
    x <- gsub(paste0("\b",w,"\b"), toupper(w), x)
  x
}
ordered_labels <- rev(c(top_up %>% arrange(NES) %>% pull(pathway) %>% clean_name(),
                         top_down %>% arrange(desc(NES)) %>% pull(pathway) %>% clean_name()))
plot_data <- plot_data %>%
  mutate(label=factor(clean_name(pathway), levels=ordered_labels),
         direction=ifelse(NES>0,"Up in APOE4","Down in APOE4"),
         sig_label=case_when(padj<0.001~"***",padj<0.01~"**",padj<0.05~"*",TRUE~""))

p <- ggplot(plot_data, aes(x=NES, y=label, fill=direction)) +
  geom_col(width=0.72, colour="white", linewidth=0.25) +
  geom_text(aes(label=sig_label, x=ifelse(NES>0,NES+0.05,NES-0.05), hjust=ifelse(NES>0,0,1)),
            size=3.8, fontface="bold", colour="grey20") +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
  geom_hline(yintercept=nrow(top_down)+0.5, linetype="dashed", colour="grey55", linewidth=0.4) +
  scale_fill_manual(values=c("Up in APOE4"="#C0392B","Down in APOE4"="#2980B9"), name=NULL) +
  scale_x_continuous(limits=c(-2.5,2.7), breaks=seq(-2,2,by=1), expand=c(0,0)) +
  labs(x="Normalized Enrichment Score (NES)", y=NULL,
       title="BIONMG APOE E4 vs E3 — Microglia-Relevant Pathway Enrichment",
       subtitle="GO:BP · Reactome · MSigDB Hallmarks  |  fgsea multilevel, padj < 0.05  |  Top 10 up + significant down") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold",size=12), plot.subtitle=element_text(size=9,colour="grey40"),
        axis.text.y=element_text(size=9.5,colour="black"), axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10), legend.position="top", legend.text=element_text(size=9.5),
        panel.grid.major.x=element_line(colour="grey90",linewidth=0.3), plot.margin=margin(10,25,10,10))

tmp <- "/workspace/BIONMG_E4vsE3_fgsea_barchart.png"
ggsave(tmp, plot=p, width=9, height=7, dpi=300, bg="white")
system(paste("cp", tmp, file.path(out_dir, "figures", "BIONMG_E4vsE3_fgsea_barchart.png")))
cat("Done.
")

