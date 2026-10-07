import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import pearsonr
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"]})
PB="../../WGCNA_10K/perdonor_module_boxplots/data"; regions=["EC","ITG","PFC","V1"]
DAM=["TREM2","TYROBP","GPNMB","CD68","CTSD","AXL","LPL","CD9","CAPG"]
rows=[]
for r in regions:
    d=pd.read_csv(f"{PB}/{r}_pseudobulk_log2cpm.csv").set_index("donor_id")
    m=pd.read_csv(f"{PB}/{r}_donor_meta.csv").set_index("donor_id")
    z=(d-d.mean())/d.std(ddof=0)
    df=pd.DataFrame(index=d.index)
    df["DAB1"]=z["DAB1"]; df["DAB2"]=z["DAB2"]; df["DAM"]=z[DAM].mean(axis=1)
    df["cerad"]=m.loc[df.index,"cerad_num"].values; rows.append(df)
A=pd.concat(rows)

fig,axes=plt.subplots(1,2,figsize=(14,6.2))
for ax,gene,title in [(axes[0],"DAB1","DAB1 expression vs DAM state"),
                      (axes[1],"DAB2","DAB2 expression vs DAM state")]:
    x=A[gene]; y=A["DAM"]; r,p=pearsonr(x,y)
    sc=ax.scatter(x,y,c=A.cerad,cmap="YlOrRd",s=95,edgecolors="#333",linewidths=0.7,zorder=3,vmin=0,vmax=3)
    b,a=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50)
    ax.plot(xs,a+b*xs,color="#2166AC",lw=2.4,zorder=2)
    ax.set_xlabel(f"{gene} expression (z-score, within region)",fontsize=12)
    ax.set_ylabel("DAM / activation signature score" if gene=="DAB1" else "",fontsize=11)
    ax.set_title(title,fontsize=13,fontweight="bold")
    ax.text(0.04,0.95,f"r = {r:+.2f}\np = {p:.1e}\nn = {len(A)}",transform=ax.transAxes,va="top",
            fontsize=12,fontweight="bold",bbox=dict(boxstyle="round,pad=0.4",fc="white",ec="#ccc"))
    ax.grid(alpha=0.25,zorder=0)
cb=fig.colorbar(sc,ax=axes,pad=0.02,fraction=0.045); cb.set_label("CERAD score",fontsize=11); cb.set_ticks([0,1,2,3])
fig.suptitle("RELN adaptors DAB1 / DAB2 vs microglial DAM state (AD Progression Atlas, per-donor pseudobulk)",
             fontsize=14,fontweight="bold",y=0.98)
out="figures/DAB1_DAB2_vs_DAM_correlation.png"; fig.savefig(out,dpi=300,facecolor="white",bbox_inches="tight")
print("saved",out)
