import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, networkx as nx, numpy as np, pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"]})
STAT={"RELN":2.69,"APOE":3.74,"APP":1.77,"MAPT":1.24,"PRNP":3.22,"CDK5R1":2.14,
 "GRIN2A":1.88,"GRIN2B":1.80,"NCSTN":2.02,"APBB1":2.22,"APBB2":2.54,"HSP90AA1":2.74,
 "SNCA":2.29,"ITGB1":1.12,"ITGB8":1.77,"LRP2":1.99,"NPC1L1":1.56,"ABCA1":1.55,
 "CUBN":1.41,"MAP1B":1.34,"PIK3CA":1.58,
 "DAB1":-0.13,"DAB2":-0.71,"VLDLR":-0.24,"LRP8":-0.53,"LRP1":-0.41,"CDK5":0.82,
 "FYN":0.06,"SORL1":0.61,"SRC":-1.14}
ALL=STAT; POS=0.50; NEG=-0.30
R=pd.read_csv("coexpr_meanR_REGIONAWARE.csv",index_col=0)
gs=[g for g in R.index if g in ALL]
pos_e=[(gs[i],gs[j],float(R.loc[gs[i],gs[j]])) for i in range(len(gs)) for j in range(i+1,len(gs))
       if pd.notna(R.loc[gs[i],gs[j]]) and R.loc[gs[i],gs[j]]>=POS]
pr=pd.read_csv("coexpr_pairs_padj_REGIONAWARE.csv")
neg_e=[(x.a,x.b,float(x.r)) for _,x in pr.iterrows() if x.a in ALL and x.b in ALL and x.r<=NEG and x.padj<0.05]
G=nx.Graph(); G.add_nodes_from(ALL)
for a,b,r in pos_e: G.add_edge(a,b)
for a,b,r in neg_e: G.add_edge(a,b)
ccs=sorted(nx.connected_components(G),key=len,reverse=True)
main=set(ccs[0]); comp=[n for n in G.nodes() if n in main]
side=[n for n in G.nodes() if n not in main]
pref=["DAB1","APBB1","CDK5","LRP8","ABCA1"]
iso_order=[n for n in pref if n in side]+sorted(x for x in side if x not in pref)
rd=cm.get_cmap("RdBu_r"); rdn=TwoSlopeNorm(0.,-4.,4.)
def lum(c): return 0.299*c[0]+0.587*c[1]+0.114*c[2]
def fe(n,s): return rd(rdn(s)),("#B2182B" if s>=0.25 else "#2166AC" if s<=-0.25 else "#AAAAAA")
Gc=G.subgraph(comp)
pc=nx.spring_layout(Gc,k=2.05,seed=4,iterations=2200)
xs=np.array([pc[n][0] for n in comp],float); ys=np.array([pc[n][1] for n in comp],float)
xs=(xs-xs.min())/(xs.max()-xs.min()+1e-9); ys=(ys-ys.min())/(ys.max()-ys.min()+1e-9)
LX0,LX1,LY0,LY1=0.11,1.47,0.13,0.90; X=LX0+xs*(LX1-LX0); Y=LY0+ys*(LY1-LY0)
MIN=0.185
for _ in range(800):
    moved=False
    for i in range(len(comp)):
        for j in range(i+1,len(comp)):
            dx=X[j]-X[i]; dy=Y[j]-Y[i]; d=(dx*dx+dy*dy)**0.5
            if d<MIN:
                if d<1e-9: dx,dy,d=0.01,0.0,0.01
                p=(MIN-d)/2.0; ux,uy=dx/d,dy/d
                X[i]-=ux*p; Y[i]-=uy*p; X[j]+=ux*p; Y[j]+=uy*p; moved=True
    X=np.clip(X,LX0,LX1); Y=np.clip(Y,LY0,LY1)
    if not moved: break
pos={n:(X[i],Y[i]) for i,n in enumerate(comp)}
SX=1.80
for k,n in enumerate(iso_order): pos[n]=(SX,0.88-k*0.155)
CM=1/2.54
fig=plt.figure(figsize=(18.72*CM,9.35*CM)); ax=fig.add_axes([0.006,0.235,0.988,0.755]); ax.set_axis_off()
ax.set_xlim(-0.02,2.00); ax.set_ylim(-0.02,1.02)
ax.add_patch(FancyBboxPatch((1.63,0.05),0.35,0.92,boxstyle="round,pad=0,rounding_size=0.015",lw=0.7,edgecolor="#C9CCD1",facecolor="#F5F6F7",zorder=0))
def ewp(r): return 0.4+1.1*(abs(r)-POS)/(0.95-POS)
def ewn(r): return 0.6+1.0*(abs(r)-0.30)/(0.5-0.30)
for a,b,r in pos_e: ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color="#AEB2B8",lw=ewp(r),alpha=0.5,zorder=1,solid_capstyle="round")
for a,b,r in neg_e: ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color="#7B4EA3",lw=ewn(r),alpha=0.85,zorder=2,ls=(0,(3,1.8)))
sideset=set(iso_order)
for a,b,r in pos_e:
    if a in sideset and b in sideset:
        ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color="#AEB2B8",lw=1.4,alpha=0.7,zorder=1,solid_capstyle="round")
FS=7.4; maxlen=max(len(n) for n in G.nodes()); Dmin=maxlen*FS*0.55+4
for n in G.nodes():
    s=ALL[n]; f,e=fe(n,s); Dd=Dmin+2.6*abs(s)
    ax.scatter([pos[n][0]],[pos[n][1]],s=Dd**2,c=[f],edgecolors=e,linewidths=1.0,zorder=3)
    ax.text(pos[n][0],pos[n][1],n,ha="center",va="center",fontsize=FS,fontweight="bold",color="white" if lum(f[:3])<0.55 else "#1a1a1a",zorder=4)
cax=fig.add_axes([0.075,0.115,0.24,0.036]); cb=fig.colorbar(ScalarMappable(norm=rdn,cmap=rd),cax=cax,orientation="horizontal")
cb.set_label("CERAD assoc. (mean Wald): blue ↓ – red ↑",fontsize=9,labelpad=2.5)
cb.set_ticks([-4,0,4]); cb.ax.tick_params(labelsize=8,length=2.5,pad=1.5)
eh=[Line2D([0],[0],color="#AEB2B8",lw=2.0,label="positive co-expression (r ≥ +0.5)"),
    Line2D([0],[0],color="#7B4EA3",lw=2.0,ls=(0,(3,1.8)),label="significant negative correlation (FDR<0.05, r ≤ −0.3)")]
ax.legend(handles=eh,loc="lower left",bbox_to_anchor=(0.34,-0.135),ncol=1,frameon=False,fontsize=9,handlelength=2.4,handletextpad=0.6,labelspacing=0.6)
fig.text(0.985,0.02,"region-aware: Fisher-z across 4 regions  ·  n = 113  ·  node size ∝ |mean Wald|",ha="right",va="bottom",fontsize=7.8,color="#8a8a8a",style="italic")
fig.savefig("Fig1d_coexpr_REGIONAWARE.png",dpi=600,facecolor="white")
print("pos",len(pos_e),"anti",len(neg_e),"side-box",iso_order)
