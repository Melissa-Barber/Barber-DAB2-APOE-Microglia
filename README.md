[README.md](https://github.com/user-attachments/files/33156558/README.md)
# Barber et al. — DAB2 / APOE microglia: analysis code

Analysis and figure-generation code for **"The adaptor protein DAB2 cell autonomously links APOE
with microglial state dysregulation relevant to Alzheimer's disease"** (Barber et al.).

Transcriptomic analysis of APOE-isogenic iPSC-microglia (DAB1/DAB2 knockdown),
human AD-brain microglia (AD Progression Atlas; Haney; Galatro; Kracht), and
weighted gene co-expression networks (WGCNA). Differential expression uses
DESeq2 / PyDESeq2; pathway enrichment uses fgsea (fgseaMultilevel) and
clusterProfiler; figures are rendered in Python (matplotlib).

Scripts are named for the figure panel they produce (all genotype contrasts are
APOE4/E4 vs APOE3/E3 throughout). Each script is standalone; helper modules
(`style.py`, `module_labels_10K.py`) sit beside the scripts that import them.

## Layout
```
_shared_DESeq2_pipeline/   counts -> PyDESeq2 -> ranked lists -> fgsea (upstream for all)
Figure1/ ... Figure7/      per-figure analysis + rendering
SuppFig1_2/  SuppFig5_6/    supplementary figures
```

## Figure -> script

| Panel | Script |
|-------|--------|
| 1a schematic | `Figure1/RELN_schematic/Fig1a_RELN_APOE_signalling_schematic.py` |
| 1b boxplots | `Figure1/Galatro_boxplots/Fig1b_Galatro_young_vs_aged_boxplots.py` |
| 1b RELN-core fgsea heatmap | `Figure1/Galatro_fgsea/Fig1b_RELNcore_fgsea_fetal_and_ageing.R` |
| 1c CERAD 4-region heatmap | `Figure1/AD_Atlas_fgsea/Fig1c_RELNcore_CERAD_4region_fgsea.R` + `Figure1/Original_CERAD_RELN_heatmap/Fig1c_RELNcore_CERAD_heatmap.py` |
| 1d co-expression network | `Figure1/RELN_pathway_CERAD_schematic/Fig1d_compute_coexpression_network.py` + `Fig1d_RELN_coexpression_network.py` |
| 1e DAB2 vs DAM scatter | `Figure1/RELN_pathway_CERAD_schematic/Fig1e_DAB2_vs_DAM_signature_scatter.py` |
| 2e genotype forest | `Figure2/Morphology_forests/Fig2e_APOE_genotype_morphology_forest.py` |
| 2f treatment forest | `Figure2/Morphology_forests/Fig2f_treatment_morphology_forest.py` |
| 3a DEGs / volcano | DEG table from `_shared_DESeq2_pipeline/.../02_run_pyDESeq2_E3E4_groupaware.py` (volcano drawn from this output) |
| 3b RELN-axis + MG-state heatmaps | `Figure3/MGstates_RELNcore_heatmap/Fig3b_E4vsE3_MGstates_RELNcore_heatmap.py` (+ `Figure3/BIONMG_CtrlsiRNA/Fig3b_E4vsE3_RELNcore_fgsea.R`, `Fig3b_E4vsE3_RELNaxis_fgsea.R`, `Fig3b_E4vsE3_RELNaxis_heatmap.py`) |
| 3c top pathways | `Figure3/BIONMG_CtrlsiRNA/Fig3c_E4vsE3_pathway_fgsea.R` |
| 4a MG-state + RELN heatmap | fgsea: `Figure4/Unbiased_HKR_v3_MULTILEVEL/Fig4_DAB2KD_pathway_fgsea.R` (rendered with the DAB-KD heatmap script below) |
| 4b pathway heatmap | `SuppFig5_6/DAB1_KD_pathway_heatmap/Fig4b_SuppFig6_DAB_KD_pathway_heatmap.py` (DAB2 panel) |
| 4c siRNA-96 vs -98 isoform dNES | `Figure7/WGCNA_10K/scripts/Fig4c_DAB2_siRNA96_vs_98_isoform_dNES.py` |
| 4d RELN-core-DAB2 vs MG-state scatter | `Figure4/RELN_core_DAB2_vs_MGstates_scatter/Fig4d_RELNcoreDAB2_vs_MGstates_scatter.py` |
| 4e blunted-response Deming scatter | drawn from the fgsea leading-edge output (no dedicated script) |
| 6 PCA (biplot + loadings) | `Figure6/NES_PCA_composite/Fig6_NES_PCA_microglial_landscape.py` |
| 7 WGCNA network build | `Figure7/WGCNA_10K/scripts/01_gene_selection.R` ... `10_module_split_5K_to_10K.R` |
| 7a module-trait heatmap | `Figure7/WGCNA_10K/scripts/Fig7a_WGCNA_module_trait_heatmap.py` |
| 7b module x MG-state overlap | `Figure7/WGCNA_10K/scripts/Fig7b_WGCNA_module_MGstate_overlap.py` |
| 7c per-module pathway ORA | `Figure7/WGCNA_10K/scripts/Fig7c_WGCNA_module_pathway_ORA.py` |
| 7e module x AD-GWAS | `Figure7/WGCNA_10K/scripts/Fig7e_WGCNA_module_AD_GWAS_enrichment.py` (+ `Fig7e_build_AD_GWAS_reference_list.py`) |
| Supp 1 RELN-axis x pathology | `SuppFig1_2/RELNaxis_genotype_pathology_correlation/SuppFig1_RELNaxis_genotype_pathology_correlation.py` |
| Supp RELN-axis forest (Haney) | `SuppFig1_2/Coray_RELNaxis_forest/SuppFig_HaneyAD_RELNaxis_forest.py` |
| Supp 4c Hammond overlap | `SuppFig1_2/Hammond_RELNcore_overlap/SuppFig4c_Hammond_RELNcore_overlap_heatmap.py` |
| Supp 6 DAB1/DAB2 KD heatmaps | `SuppFig5_6/DAB1_KD_pathway_heatmap/Fig4b_SuppFig6_DAB_KD_pathway_heatmap.py` |
| Supp 6c DAB1 pooled + complement | `SuppFig5_6/DAB1_KD_pathway_heatmap/SuppFig6c_DAB1_pooled_complement_heatmap.py` |
| Supp 6 KEGG | `SuppFig5_6/DAB1_KD_pathway_heatmap/SuppFig6_DAB1_KD_KEGG_fgsea.R` |
| Supp 5 volcanoes | `SuppFig5_6/DAB_KD_volcanoes/SuppFig5_DAB_KD_volcanoes.R` |
| Supp 5 isoform concordance | `SuppFig5_6/DAB2_isoform_genotype_scatter/SuppFig5_DAB2_isoform_concordance.R` |

Note: `Fig4b_SuppFig6_DAB_KD_pathway_heatmap.py` renders both the Fig 4b (DAB2 KD)
and Supplementary Fig 6 (DAB1 KD) pathway heatmaps. `Fig4c_...` lives under the
Figure7 folder because it shares that folder's helper modules.

**No code (standard software):** Figs 2a-d (qRT-PCR, immunostaining, imaging PCA),
Fig 5 (Reln-KO mouse histology), Fig 7d TREM2 immunostaining, and the Fig 3a/4e
volcano & Deming scatters drawn directly from the DESeq2/fgsea output tables.

## Paths - edit before running
Set absolute paths to your local copy: placeholders `<EDIT_PROJECT_ROOT>`,
`<EDIT_INPUT_DIR>`, `<EDIT_OUTPUT_DIR>`, and the `PROJECT`/path variable near the
top of some R scripts (`/Users/...` or `/Volumes/...`). Helper modules
`style.py` and `module_labels_10K.py` are bundled beside the scripts that import
them. Input data accompany the manuscript as Source/Supplementary Data; raw
sequencing data are at GEO ([accession]).

## Software
- **R >= 4.3**: DESeq2, fgsea, msigdbr, WGCNA, clusterProfiler, data.table, dplyr, readr
- **Python >= 3.10**: pandas, numpy, matplotlib, pydeseq2, gseapy, scipy, statsmodels, networkx, adjustText

## Citation
Please cite Barber et al. and the archived release (Zenodo DOI: 10.5281/zenodo.XXXXXX).

## License
MIT - see LICENSE.
