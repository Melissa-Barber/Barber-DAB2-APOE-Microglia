import re,pandas as pd,numpy as np
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":["Arial","DejaVu Sans"]})
D="<EDIT_INPUT_DIR>/FINAL FIGURES & SCRIPTS_JUNE2026/Figures/Claude Python scripts/Galatro_boxplots/data"
META=f"{D}/GSE99074-GPL16791_series_matrix.txt"; VOOM=f"{D}/GSE99074_HumanMicrogliaBrainVoomNormalization.txt"
# CORRECTED ids (DAB2 fixed)
ENS={"LRP8":"ENSG00000157193","VLDLR":"ENSG00000147852","LRP1":"ENSG00000123384","LRP2":"ENSG00000081479",
 "ITGB1":"ENSG00000150093","ITGB8":"ENSG00000105855","RELN":"ENSG00000189056","APOE":"ENSG00000130203",
 "DAB1":"ENSG00000173406","DAB2":"ENSG00000153071"}
GROUPS={"Adults":(34,60),"Aged":(60,103)}
def grp(a):
    for g,(lo,hi) in GROUPS.items():
        if lo<=a<hi: return g
    return None
titles,ages=[],[]
for line in open(META):
    if line.startswith("!Sample_title"): titles=[t.strip().strip('"') for t in line.strip().split("\t")[1:]]
    elif line.startswith("!Sample_characteristics_ch1"):
        vals=[t.strip().strip('"') for t in line.strip().split("\t")[1:]]
        if vals and vals[0].lower().startswith("age:"):
            ages=[int(re.search(r"(\d+)",v).group(1)) if re.search(r"(\d+)",v) else None for v in vals]
meta=pd.DataFrame({"sample_id":titles,"age":ages}).dropna(subset=["age"]); meta["age"]=meta.age.astype(int)
meta["g"]=meta.age.apply(grp); meta=meta.dropna(subset=["g"])
voom=pd.read_csv(VOOM,sep="\t",index_col=0)
keep=[s for s in meta.sample_id if s in voom.columns]; meta=meta[meta.sample_id.isin(keep)]
young=meta[meta.g=="Adults"].sample_id.tolist(); aged=meta[meta.g=="Aged"].sample_id.tolist()
# stats all 10 genes (corrected DAB2)
rows=[]
for g,e in ENS.items():
    if e not in voom.index: rows.append(dict(gene=g,MW_p=np.nan)); continue
    y=voom.loc[e,young].astype(float).values; a=voom.loc[e,aged].astype(float).values
    frac=float((np.concatenate([y,a])>0).mean())
    rows.append(dict(gene=g,n_young=len(y),n_aged=len(a),median_young=round(np.median(y),3),
                     median_aged=round(np.median(a),3),MW_p=mannwhitneyu(y,a,alternative="two-sided").pvalue,
                     detection_frac=round(frac,3),filter_status="pass" if frac>=0.8 else "fail"))
st=pd.DataFrame(rows)
ok=st.filter_status=="pass"; p=st.loc[ok,"MW_p"].values; o=p.argsort(); r=np.empty(len(p)); r[o]=np.arange(1,len(p)+1)
q=np.minimum.accumulate((p*len(p)/r)[o][::-1])[::-1]; qb=np.empty(len(p)); qb[o]=np.clip(q,0,1)
st.loc[ok,"MW_q_BH"]=qb
st.to_csv("<EDIT_OUTPUT_DIR>/Galatro_10genes_MW_stats_DAB2corrected.csv",index=False)
print(st[["gene","median_young","median_aged","MW_p","MW_q_BH"]].to_string(index=False))

# render Fig 1b: 6 genes, 2x3
PLOT=[["LRP8","VLDLR","ITGB1"],["LRP1","DAB1","DAB2"]]
COL={"Adults":"#2E6FB0","Aged":"#D7263D"}
qmap=dict(zip(st.gene,st.MW_q_BH)); pmap=dict(zip(st.gene,st.MW_p))
def star(q): return "ns" if (pd.isna(q) or q>=0.05) else ("***" if q<0.001 else "**" if q<0.01 else "*")
fig,axes=plt.subplots(2,3,figsize=(6.6,4.6))
np.random.seed(1)
yl=int(meta[meta.g=="Adults"].age.min()); yh=int(meta[meta.g=="Adults"].age.max())
al=int(meta[meta.g=="Aged"].age.min()); ah=int(meta[meta.g=="Aged"].age.max())
xlab=[f"Young\n({yl}-{yh} y)\nn={len(young)}",f"Aged\n({al}-{ah} y)\nn={len(aged)}"]
for ri,row in enumerate(PLOT):
    for ci,g in enumerate(row):
        ax=axes[ri,ci]; e=ENS[g]
        data=[voom.loc[e,young].astype(float).values,voom.loc[e,aged].astype(float).values]
        for i,(grpname,v) in enumerate(zip(["Adults","Aged"],data)):
            ax.boxplot([v],positions=[i],widths=0.55,showfliers=False,patch_artist=True,
                boxprops=dict(facecolor=COL[grpname],alpha=.45,edgecolor=COL[grpname]),
                medianprops=dict(color="k"),whiskerprops=dict(color=COL[grpname]),capprops=dict(color=COL[grpname]))
            ax.scatter(np.random.normal(i,.07,len(v)),v,s=7,c=COL[grpname],alpha=.8,edgecolor="white",lw=.3,zorder=3)
        allv=np.concatenate(data); yspan=allv.max()-allv.min(); top=allv.max()+0.12*yspan
        ax.plot([0,0,1,1],[top-0.03*yspan,top,top,top-0.03*yspan],lw=0.8,c="#555")
        ax.text(0.5,top+0.01*yspan,star(qmap.get(g)),ha="center",va="bottom",fontsize=8,color="#555")
        ax.set_title(g,fontsize=10,fontweight="bold",style="italic",pad=3)
        ax.set_xticks([0,1]); ax.set_xticklabels(xlab,fontsize=6.5)
        ax.set_ylim(top=top+0.16*yspan); ax.tick_params(labelsize=7.5)
        if ci==0: ax.set_ylabel("voom log$_2$ CPM",fontsize=8)
        ax.spines[["top","right"]].set_visible(False)
fig.suptitle("Fig 1b (corrected: DAB2 = ENSG00000153071)  —  Galatro young vs aged microglia",fontsize=9,y=1.0)
plt.tight_layout()
plt.savefig("<EDIT_OUTPUT_DIR>/Fig1b_Galatro_young_vs_aged_DAB2corrected.png",dpi=300,bbox_inches="tight",facecolor="white")
print("saved fig")
