import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"]})
B="<EDIT_INPUT_DIR>/FINAL FIGURES & SCRIPTS_JUNE2026/Figures/Claude revised scripts/RELN_core_DAB12_vs_MGstates_scatter/data"
master=pd.read_csv(f"{B}/NES_padj_RELN_and_MGstates_all_contrasts.csv")
reln=pd.read_csv(f"{B}/RELN_core_DAB12_fgsea_per_contrast.csv")
AD={"AD Progression Atlas","Coray AD snRNA","iPSC-MG APOE","BIONMG genotype","BIONMG KD"}
mg=master[master.contrast_group.isin(AD)]
piv=mg.pivot_table(index=["contrast_id","contrast_group"],columns="pathway",values="NES",aggfunc="first").reset_index()
rp=reln.pivot_table(index="contrast_id",columns="pathway",values="NES",aggfunc="first").reset_index()
m=piv.merge(rp,on="contrast_id",how="left")
m=m[~m.contrast_id.str.contains("_E2")].copy()   # REMOVE APOE2  -> n=20
X="RELN_core_DAB2"
# colour by data SOURCE; Haney human brain = GREEN, Haney iPSC-APOE = YELLOW
GCOL={"AD Progression Atlas":"#D7263D","Coray AD snRNA":"#2E8B57","iPSC-MG APOE":"#F5B041",
      "BIONMG genotype":"#7FB3D5","BIONMG KD":"#1F4E79"}
GLAB={"AD Progression Atlas":"AD Progression Atlas (Wachter 2024)",
      "Coray AD snRNA":"Haney AD human brain, PFC (2024)",
      "iPSC-MG APOE":"Haney APOE iPSC-MG (2024)",
      "BIONMG genotype":"iPSC APOE genotype, Ctrl-siRNA (this study)",
      "BIONMG KD":"iPSC APOE DAB1/DAB2-KD (this study)"}
GORDER=["AD Progression Atlas","Coray AD snRNA","iPSC-MG APOE","BIONMG genotype","BIONMG KD"]
MG=[("SalaFrigerio2019_ARM","ARM\n(SalaFrigerio)"),("Homeostatic (HM)","Homeostatic\n(consensus)"),
    ("Mathys2019_MG_Homeostatic","Homeostatic\n(Mathys)"),("Haney2024_LDAM_markers","LDAM\n(Haney)"),
    ("Cytokines response 1 (CRM-1)","CRM-1"),("Cytokines response 2 (CRM-2)","CRM-2")]
E4E3={"BIONMG_E4_vs_E3_Ctrl","Coray_E44_vs_E33_LDneg","Coray_E44_vs_E33_LDpos","Coray_snMG_AD44_vs_AD33","MIC_E4_vs_E3"}
def bh(p):
    p=np.array(p); o=p.argsort(); r=np.empty(len(p)); r[o]=np.arange(1,len(p)+1)
    q=p*len(p)/r; qs=np.minimum.accumulate(q[o][::-1])[::-1]; out=np.empty(len(p)); out[o]=np.clip(qs,0,1); return out
ps=[spearmanr(m[[X,c]].dropna()[X],m[[X,c]].dropna()[c])[1] for c,_ in MG]; fdr=bh(ps)
fig,axes=plt.subplots(1,6,figsize=(27,6.8),sharex=True,sharey=True)
xall=m[X].dropna(); xmin,xmax=xall.min()-0.3,xall.max()+0.45
for ax,(col,lab),q in zip(axes,MG,fdr):
    sub=m[[X,col,"contrast_group","contrast_id"]].dropna().copy()
    x=sub[X].values; y=sub[col].values; rho,_=spearmanr(x,y)
    mm,bb=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50); ax.plot(xs,mm*xs+bb,color="#333",lw=1.4,zorder=2)
    rng=np.random.RandomState(0); boots=[]
    for _ in range(400):
        ii=rng.randint(0,len(x),len(x))
        try: s2,i2=np.polyfit(x[ii],y[ii],1); boots.append(s2*xs+i2)
        except: pass
    if boots:
        ba=np.array(boots); lo,hi=np.percentile(ba,[2.5,97.5],axis=0); ax.fill_between(xs,lo,hi,color="#999",alpha=0.15,zorder=1)
    for grp in GORDER:
        c=GCOL[grp]
        for _,r in sub[sub.contrast_group==grp].iterrows():
            mk="^" if r.contrast_id in E4E3 else "o"
            ax.scatter(r[X],r[col],s=120 if mk=="^" else 80,c=c,marker=mk,edgecolor="#111",linewidth=0.7,alpha=0.95,zorder=3)
    sym="**" if q<0.01 else ("*" if q<0.05 else "")
    ax.text(0.04,0.97,f"ρ={rho:+.2f}{sym}\nFDR={q:.3f}\nn={len(sub)}",transform=ax.transAxes,ha="left",va="top",
            fontsize=12,fontweight="bold",bbox=dict(fc="white",ec="none",alpha=0.85,pad=2))
    ax.set_title(lab,fontsize=13,fontweight="bold",pad=6)
    ax.axhline(0,color="#bbb",lw=0.6,ls="--"); ax.axvline(0,color="#bbb",lw=0.6,ls="--")
    ax.set_xlim(xmin,xmax); ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("MG state NES",fontsize=14,fontweight="bold")
fig.text(0.45,0.02,"RELN-core-DAB2 (113g) NES",ha="center",fontsize=14,fontweight="bold")
fig.suptitle("Fig 4d — RELN-core-DAB2 vs AD-relevant microglial states  (APOE2 removed, n=20)",fontsize=14,fontweight="bold",y=1.01)
handles=[Line2D([0],[0],marker="o",color="w",markerfacecolor=GCOL[g],markeredgecolor="#111",markersize=12,label=GLAB[g]) for g in GORDER]
handles+=[Line2D([0],[0],marker="^",color="w",markerfacecolor="#BBB",markeredgecolor="#111",markersize=13,label="APOE4/E4 vs APOE3/E3 (triangles)"),
          Line2D([0],[0],marker="o",color="w",markerfacecolor="#BBB",markeredgecolor="#111",markersize=12,label="all other contrasts (circles)")]
fig.legend(handles=handles,loc="center left",bbox_to_anchor=(0.905,0.5),fontsize=10.5,frameon=False,title="Contrast group / marker",title_fontsize=11)
fig.subplots_adjust(right=0.90,bottom=0.12,top=0.84,wspace=0.12)
out="<EDIT_OUTPUT_DIR>/Fig4d_RELNcoreDAB2_vs_MGstates_noAPOE2_Haneygreen.png"
fig.savefig(out,dpi=300,bbox_inches="tight",facecolor="white"); plt.close(fig)
print("saved",out)
# report which contrasts are green (Haney human brain) and yellow (Haney iPSC)
print("GREEN (Haney human brain):",sorted(m[m.contrast_group=="Coray AD snRNA"].contrast_id.unique()))
print("YELLOW (Haney iPSC):",sorted(m[m.contrast_group=="iPSC-MG APOE"].contrast_id.unique()))
