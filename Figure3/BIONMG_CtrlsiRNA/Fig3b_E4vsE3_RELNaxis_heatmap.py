#!/usr/bin/env python3
# Final E4 vs e3 heatmap (WALD-ranked), approved design: all-black text, lightened
# cells, fat two-line cells (NES/stars on top, leading-edge genes below), 10.77 cm.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd, numpy as np
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"]})

CSV="<EDIT_INPUT_DIR>/fgsea_E4vsE3_axis_heatmap_WALD/data/E4vsE3_axis_heatmap_table.csv"
df=pd.read_csv(CSV)

def clean(p):
    import re
    p=re.sub(r"^(GOBP_|REACTOME_|HALLMARK_|KEGG_(LEGACY_)?)","",p)
    s=p.replace("_"," ").title()
    small={"Of","The","Via","In","To","And","For","By","On","With"}
    s=" ".join(w if (i==0 or w not in small) else w.lower() for i,w in enumerate(s.split()))
    for a,b in [("Ecm","ECM"),("Tnfa","TNFα"),("Nfkb","NF-κB"),("Dna","DNA"),
                ("Apoe","APOE"),("Tlr","TLR"),("Mhc","MHC")]:
        s=re.sub(r"\b"+a+r"\b",b,s)
    return s
SHORT={
 "Cell Morphogenesis Involved in Neuron Differentiation":"Cell Morphogenesis in Neuron Differentiation",
 "Positive Regulation of Ubiquitin Dependent Protein Catabolic Process":"Pos. Reg. Ubiquitin-Dependent Catabolism",
 "Antigen Receptor Mediated Signaling Pathway":"Antigen Receptor-Mediated Signaling",
 "Immune Response Regulating Signaling Pathway":"Immune Response-Regulating Signaling",
 "Toll Like Receptor Signaling Pathway":"Toll-Like Receptor Signaling",
 "Autophagosome Lysosome Fusion":"Autophagosome–Lysosome Fusion",
}
df["name"]=[SHORT.get(clean(p),clean(p)) for p in df["pathway"]]
df["gene_lab"]=df["gene_lab"].fillna("—")

up=df[df.NES>0].sort_values("NES")           # ascending -> smallest at top
dn=df[df.NES<0].sort_values("NES",ascending=False)  # -1.46 .. -2.25
rows=list(up.itertuples())+list(dn.itertuples())
n_up=len(up)

# lightened RdBu so black text reads
_b=cm.get_cmap("RdBu_r"); _c=_b(np.linspace(0,1,256)); _c[:,:3]=_c[:,:3]*0.58+0.42
rd=LinearSegmentedColormap.from_list("l",_c); norm=TwoSlopeNorm(0.,-3.,3.); INK="#1a1a1a"

fig=plt.figure(figsize=(19/2.54,10.77/2.54))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_axis_off()

ax.text(0.5,0.966,"BIONMG iPSC-MG Ctrl-siRNA  —  APOE-ε4 vs ε3 pathway enrichment",
        ha="center",va="center",fontsize=9.6,fontweight="bold",color=INK)
ax.text(0.5,0.934,"Hallmark · KEGG · Reactome · GO:BP  |  multilevel fgsea (min 15, max 500)  ·  DESeq2 Wald-stat ranking",
        ha="center",va="center",fontsize=5.9,style="italic",color="#6b7078")
ax.text(0.5,0.908,"italic = RELN/APOE-axis genes in leading edge (k of N = set size)   ·   *** <0.001   ** <0.01   * <0.05",
        ha="center",va="center",fontsize=5.9,color="#6b7078")

TOPy,BOTy=0.862,0.040
N=len(rows); rh=(TOPy-BOTy)/N
def yc(i): return TOPy-rh*(i+0.5)
NAMEx=0.360; CX0,CX1=0.374,0.740

for i,r in enumerate(rows):
    y=yc(i); f=rd(norm(r.NES))
    ax.add_patch(FancyBboxPatch((CX0,y-rh*0.44),CX1-CX0,rh*0.88,
        boxstyle="round,pad=0,rounding_size=0.006",lw=0.7,edgecolor="white",facecolor=f,zorder=2))
    ax.text((CX0+CX1)/2,y+rh*0.21,f"NES = {r.NES:+.2f}    {r.stars}",ha="center",va="center",
            fontsize=6.4,fontweight="bold",color=INK,zorder=3)
    ax.text((CX0+CX1)/2,y-rh*0.21,f"{r.gene_lab}   ({r.k_axis}/{r.size})",ha="center",va="center",
            fontsize=5.8,style="italic",fontweight="bold",color=INK,zorder=3)
    ax.text(NAMEx,y,r.name,ha="right",va="center",fontsize=6.2,color=INK,zorder=3)

ysep=(yc(n_up-1)+yc(n_up))/2
ax.plot([0.02,CX1],[ysep,ysep],color="#9aa0a6",lw=0.7,ls=(0,(4,2)),zorder=1)
ax.text(0.758,yc((n_up-1)/2),"Up in APOE4",ha="center",va="center",rotation=90,fontsize=6.2,color="#B2182B",fontweight="bold")
ax.text(0.758,yc(n_up+(N-n_up-1)/2),"Down",ha="center",va="center",rotation=90,fontsize=6.2,color="#2166AC",fontweight="bold")

cax=fig.add_axes([0.82,0.12,0.013,0.58])
cb=fig.colorbar(ScalarMappable(norm=norm,cmap=rd),cax=cax)
cb.set_label("NES",fontsize=6.5,labelpad=1.5); cb.set_ticks([-3,0,3]); cb.ax.tick_params(labelsize=5.8,length=2,pad=1)

fig.savefig("E4vsE3_axis_heatmap_WALD.png",dpi=600,facecolor="white")
print("saved; rows:",N,"up:",n_up)
