#!/usr/bin/env python3
"""Genotype (E4-E3) morphology/marker forest — significant hits (FDR<0.05).
Rows grouped by readout (Length, Area, Roundness, then markers) with treatments
ordered Ctrl -> LPS -> Myelin -> IFNg. Ctrl rows are coloured black; other
conditions share one colour. No FDR legend (stars mark significance).
x-axis: 'E4 - E3 (log10)'. Reads ../data/genotype_forest_results.csv."""
import os, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
CSV=os.path.join(ROOT,"data","genotype_forest_results.csv")
OUT=os.path.join(ROOT,"figures","genotype_forest_significant.png"); os.makedirs(os.path.dirname(OUT),exist_ok=True)
plt.rcParams.update({"font.size":9,"font.family":["Arial","Liberation Sans","Helvetica","DejaVu Sans"],
    "axes.titlesize":9,"axes.labelsize":9,"xtick.labelsize":9,"ytick.labelsize":9,
    "axes.spines.top":False,"axes.spines.right":False,"axes.grid":False})
READOUT_ORDER=["Length","Area","Roundness","DAB1","DAB1pY232","DAB2","CAPG","BODIPY"]
TREAT_ORDER=["Ctrl","LPS","Myelin","IFNg"]
BASE_COL="#15406b"; CTRL_COL="#000000"
res=pd.read_csv(CSV); sig=res[res.fdr<0.05].copy()
sig["_rk"]=sig.readout.map(lambda r: READOUT_ORDER.index(r) if r in READOUT_ORDER else 99)
sig["_tk"]=sig.treatment.map(lambda t: TREAT_ORDER.index(t) if t in TREAT_ORDER else 99)
sig=sig.sort_values(["_rk","_tk"]).reset_index(drop=True)  # top -> bottom
n=len(sig)
fig,ax=plt.subplots(figsize=(5.6,5.97/2.54))
for i,(_,x) in enumerate(sig.iterrows()):
    y=n-1-i; c=CTRL_COL if x.treatment=="Ctrl" else BASE_COL
    ax.errorbar(x.est,y,xerr=[[x.est-x.lo],[x.hi-x.est]],fmt='o',ms=4.5,color=c,ecolor=c,elinewidth=1.2,capsize=2.5,zorder=3)
    st="***" if x.fdr<.001 else "**" if x.fdr<.01 else "*"
    ax.text(x.hi+.006,y,st,va="center",fontsize=9,color=c)
ax.axvline(0,ls="--",color="k",lw=.6)
ax.set_yticks([n-1-i for i in range(n)]); ax.set_yticklabels([f"{x.readout} · {x.treatment}" for _,x in sig.iterrows()])
ax.set_ylim(-0.6,n-0.4)
ax.set_xlabel("E4 − E3  (log10)")
ax.set_title("Significant APOE4 genotype effects")
plt.tight_layout(pad=0.4); plt.savefig(OUT,dpi=300); print("saved",OUT)
