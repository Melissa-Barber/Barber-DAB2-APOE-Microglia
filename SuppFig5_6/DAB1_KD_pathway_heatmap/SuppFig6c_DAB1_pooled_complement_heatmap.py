#!/usr/bin/env python3
"""Supp 6c — DAB1 KD pathway heatmap switched to POOLED siRNA-893+894 genotype
columns, with REACTOME_COMPLEMENT_CASCADE added (significant: APOE4-pooled **).
Replicates render_BIONMG_DAB1_DAB2_KD_heatmap.py render() + _shared.style exactly."""
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# ---- _shared.style (verbatim) ----
CM = 1/2.54; DPI = 600
FONTS = {"title":14,"value":20,"star":18,"row":18,"col":14,"caption":9,"cbar_lbl":13,"cbar_tick":11}
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","Liberation Sans","DejaVu Sans"],
    "font.size":10.0,"axes.linewidth":0.6,"savefig.dpi":DPI,"savefig.bbox":"tight","pdf.fonttype":42,"ps.fonttype":42})
NES_CMAP = LinearSegmentedColormap.from_list("div_bwr", ["#2166AC","#FFFFFF","#B2182B"], N=256)
NES_NORM = TwoSlopeNorm(vcenter=0.0, vmin=-2, vmax=2)
def stars(p):
    if p is None or (isinstance(p,float) and p!=p): return ""
    return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""

# ---- inputs ----
U="<EDIT_INPUT_DIR>"
CSV=f"{U}/Figures/BIONMG_DAB1siRNA_KD/data/BIONMG_DAB1_KD_3contrasts_fgsea_full.csv"
RELN=set()
for f in [f"{U}/Final Figures  June62026/Figures/gene_sets/RELN_core_DAB1_115genes.csv",
          f"{U}/Final Figures  June62026/Figures/gene_sets/RELN_core_DAB2_113genes.csv"]:
    d=pd.read_csv(f); col="gene_symbol" if "gene_symbol" in d.columns else d.columns[0]
    RELN.update(d[col].dropna().astype(str).str.upper().str.strip())
print("RELN query genes:",len(RELN))

CONTRASTS=["E4_pooled_DAB1","E3_pooled_DAB1","int_893"]
PID2LABEL={
 "REACTOME_CHOLESTEROL_BIOSYNTHESIS":"Cholesterol Biosynthesis",
 "HALLMARK_FATTY_ACID_METABOLISM":"Fatty Acid Metabolism",
 "HALLMARK_OXIDATIVE_PHOSPHORYLATION":"Oxidative Phosphorylation",
 "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION":"Mitochondrial Proteostasis",
 "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION":"ECM Organization",
 "HALLMARK_MTORC1_SIGNALING":"mTORC1 Signaling",
 "GOBP_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS":"Ubiquitin-Dependent Proteolysis",
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB":"TNFα Signaling via NF-κB",
 "HALLMARK_INFLAMMATORY_RESPONSE":"Inflammatory Response",
 "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP":"SASP",
 "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE":"DNA Damage-Induced Senescence",
 "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY":"Cytokine-Mediated Signaling",
 "GOBP_AMYLOID_BETA_CLEARANCE":"Amyloid-β Clearance",
 "GOBP_PHAGOCYTOSIS":"Phagocytosis",
 "GOBP_PHAGOCYTOSIS_ENGULFMENT":"Phagocytosis Engulfment",
 "GOBP_LYSOSOMAL_TRANSPORT":"Lysosomal Transport",
 "REACTOME_COMPLEMENT_CASCADE":"Complement Cascade",   # <-- added
}
df=pd.read_csv(CSV)
DATA={}
for pid,label in PID2LABEL.items():
    cells=[]
    for c in CONTRASTS:
        r=df[(df.pathway==pid)&(df.contrast==c)]
        if r.empty: cells.append((np.nan,"",0,0)); continue
        r=r.iloc[0]
        le={g.strip().upper() for g in str(r["leadingEdge"]).split(";") if g.strip()}
        cells.append((float(r["NES"]), stars(float(r["padj"])), len(le&RELN), int(r["size"])))
    DATA[label]=cells

# ---- ROWS (complement added to the phagocytic/endolysosomal 'manual' block) ----
ROWS=[("Cholesterol Biosynthesis","interact_up"),("Fatty Acid Metabolism","interact_up"),
 ("Oxidative Phosphorylation","interact_up"),("Mitochondrial Proteostasis","interact_up"),
 ("ECM Organization","interact_up"),("mTORC1 Signaling","interact_up"),
 ("Ubiquitin-Dependent Proteolysis","interact_dn"),("TNFα Signaling via NF-κB","interact_dn"),
 ("Inflammatory Response","interact_dn"),("SASP","interact_dn"),
 ("DNA Damage-Induced Senescence","interact_dn"),("Cytokine-Mediated Signaling","interact_dn"),
 ("Amyloid-β Clearance","manual"),("Phagocytosis","manual"),
 ("Phagocytosis Engulfment","manual"),("Complement Cascade","manual"),("Lysosomal Transport","manual")]
COL_LABELS=["APOE4 pooled\nDAB1 KD\n(siRNA-893+894)","APOE3 pooled\nDAB1 KD\n(siRNA-893+894)",
            "APOE4 vs APOE3\nDAB1 interact\n(siRNA-893)"]
n=len(ROWS)
BLOCK_GAP=0.45
def y_for(i):
    b=ROWS[i][1]
    return i if b=="interact_up" else i+BLOCK_GAP if b=="interact_dn" else i+2*BLOCK_GAP
y_max=y_for(n-1)+0.6
CELL_W,CELL_H=3.2,1.5
fig_w_cm=26.0; fig_h_cm=9.0+1.10*y_max
fig=plt.figure(figsize=(fig_w_cm*CM,fig_h_cm*CM))
ax=fig.add_axes([0.27,0.05,0.50,0.72])
for i,(pid,_b) in enumerate(ROWS):
    yc=y_for(i)*CELL_H
    for j,(v,star,nr,gs) in enumerate(DATA[pid]):
        xc=j*CELL_W
        face=NES_CMAP(NES_NORM(np.clip(v,-3,3))) if not np.isnan(v) else (0.93,0.93,0.93,1)
        ax.add_patch(Rectangle((xc-CELL_W/2,yc-CELL_H/2),CELL_W,CELL_H,facecolor=face,edgecolor="white",linewidth=1.2))
        tcol="white" if (not np.isnan(v) and abs(v)>1.85) else "#222"
        ax.text(xc,yc-0.20,f"{v:+.2f}" if not np.isnan(v) else "—",ha="center",va="center",fontsize=13,color=tcol,weight=("bold" if star else "normal"))
        if star: ax.text(xc,yc+0.40,star,ha="center",va="center",fontsize=11,color=tcol,weight="bold")
# dividers (block boundaries: up|dn after idx5, dn|manual after idx11)
for a,b,ls,cc,lw in [(5,6,"--","#666",0.8),(11,12,"-","#222",0.9)]:
    yy=(y_for(a)+y_for(b))/2*CELL_H
    ax.plot([-CELL_W/2-0.1,(len(COL_LABELS)-1)*CELL_W+CELL_W/2+0.1],[yy,yy],linestyle=ls,color=cc,linewidth=lw)
yticks=[y_for(i)*CELL_H for i in range(n)]
ylabs=[f"{ROWS[i][0]} ({DATA[ROWS[i][0]][0][2]}/{DATA[ROWS[i][0]][0][3]})" for i in range(n)]
ax.set_yticks(yticks); ax.set_yticklabels(ylabs,fontsize=FONTS["row"])
ax.set_xticks([j*CELL_W for j in range(len(COL_LABELS))]); ax.set_xticklabels(COL_LABELS,fontsize=FONTS["col"],weight="bold")
ax.xaxis.set_label_position("top"); ax.xaxis.tick_top()
ax.set_ylim(y_max*CELL_H+0.2,-CELL_H/2-0.2)
ax.tick_params(top=False,bottom=False,left=False,right=False,pad=14)
for sp in ax.spines.values(): sp.set_visible(False)
bracket_x=(len(COL_LABELS)-1)*CELL_W+CELL_W/2+0.35; text_x=bracket_x+0.20
def yrange(idxs): return min(y_for(i) for i in idxs)*CELL_H-CELL_H/2+0.05, max(y_for(i) for i in idxs)*CELL_H+CELL_H/2-0.05
groups=[([0,1,2,3,4,5],"Lipid, Mitochondria,\nECM & mTOR\n(interact ↑)","#2C2C2C"),
        ([6,7,8,9,10,11],"Inflammation &\nSenescence\n(interact ↓)","#2C2C2C"),
        ([12,13,14,15,16],"Amyloid-β, Phagocytosis,\nComplement & Endolysosomal","#555555")]
ax.set_xlim(-CELL_W/2-0.2,(len(COL_LABELS)-1)*CELL_W+CELL_W/2+2.4)
for idxs,lbl,col in groups:
    y0,y1=yrange(idxs); ax.plot([bracket_x,bracket_x],[y0,y1],color=col,linewidth=1.4)
    ax.text(text_x,(y0+y1)/2,lbl,ha="left",va="center",fontsize=FONTS["row"]-1,weight="bold",color=col,linespacing=1.15,clip_on=False)
cax=fig.add_axes([0.91,0.15,0.018,0.50])
cb=mpl.colorbar.ColorbarBase(cax,cmap=NES_CMAP,norm=NES_NORM,ticks=[-3,-2,-1,0,1,2,3])
cb.set_label("NES",fontsize=FONTS["cbar_lbl"],weight="bold"); cb.ax.tick_params(labelsize=FONTS["cbar_tick"]); cb.outline.set_linewidth(0.4)
fig.text(0.5,0.965,"BIONMG iPSC-MG — DAB1 KD Pathway Enrichment (strict-filter v3 ranking)",ha="center",va="top",fontsize=FONTS["title"]+2,weight="bold")
fig.text(0.5,0.915,"APOE4 / APOE3 pooled DAB1 knockdown (siRNA-893+894) and the genotype × DAB1-KD interaction (siRNA-893).\n"
    "fgsea multilevel; DESeq2 with ≥10 counts in ≥80% Ctrl prevalence filter (v3 ranking)\n"
    "* padj<0.05  ** <0.01  *** <0.001.   Parentheses: RELN-core-DAB1/2 leading-edge / pathway size.",
    ha="center",va="top",fontsize=FONTS["caption"],color="#555",linespacing=1.4)
out="<EDIT_OUTPUT_DIR>/SuppFig6c_DAB1_pooled_pathway_heatmap_withComplement.png"
fig.savefig(out,dpi=DPI,bbox_inches="tight",facecolor="white"); plt.close(fig)
print("saved",out)
print("complement row:",DATA["Complement Cascade"])
