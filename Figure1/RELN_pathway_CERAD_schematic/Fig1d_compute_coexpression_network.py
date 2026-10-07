import pandas as pd, numpy as np, itertools, math
from math import erf, sqrt, atanh, tanh
U="<EDIT_INPUT_DIR>/Revised Figures_August2026/RELN_pathway_CERAD_schematic/_recompute_tmp"
regs={r:pd.read_csv(f"{U}/{r}_sub.csv").set_index("donor_id") for r in ["EC","ITG","PFC","V1"]}
genes=list(regs["EC"].columns)
def phi(z): return 0.5*(1+erf(z/sqrt(2)))
def pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); 
    if a.std()==0 or b.std()==0: return np.nan
    return np.corrcoef(a,b)[0,1]

rows=[]
allcat={g: pd.concat([regs[r][g] for r in regs]).values for g in genes}  # pooled 113
for a,b in itertools.combinations(genes,2):
    zk=[]; wk=[]; rk=[]
    for r in regs:
        x=regs[r][a].values; y=regs[r][b].values; n=len(x)
        rr=pearson(x,y)
        if np.isnan(rr): continue
        rr=min(max(rr,-0.999),0.999)
        zk.append(atanh(rr)); wk.append(n-3); rk.append(rr)
    W=sum(wk); zbar=sum(w*z for w,z in zip(wk,zk))/W
    se=1/sqrt(W); Z=zbar/se; p_ra=2*(1-phi(abs(Z))); r_ra=tanh(zbar)
    # naive pooled
    rp=pearson(allcat[a],allcat[b]); npool=len(allcat[a])
    tp=rp*sqrt((npool-2)/max(1e-12,1-rp*rp)); 
    # pooled p via normal approx on fisher z
    zp=atanh(min(max(rp,-0.999),0.999)); Zp=zp*sqrt(npool-3); p_pool=2*(1-phi(abs(Zp)))
    rows.append(dict(a=a,b=b,r_regionaware=r_ra,p_regionaware=p_ra,
                     r_pooled=rp,p_pooled=p_pool,r_meanSimple=np.mean(rk)))
df=pd.DataFrame(rows)
# BH FDR on region-aware
def bh(p):
    p=np.asarray(p); n=len(p); order=np.argsort(p); ranks=np.empty(n,int); ranks[order]=np.arange(1,n+1)
    q=p*n/ranks; # enforce monotonicity
    qs=q[order]; qs=np.minimum.accumulate(qs[::-1])[::-1]; out=np.empty(n); out[order]=qs
    return np.minimum(out,1)
df["padj_regionaware"]=bh(df["p_regionaware"].values)
df["padj_pooled"]=bh(df["p_pooled"].values)
df.to_csv("coexpr_pairs_REGIONAWARE.csv",index=False)

# region-aware mean-R matrix
R=pd.DataFrame(np.eye(len(genes)),index=genes,columns=genes)
for _,x in df.iterrows():
    R.loc[x.a,x.b]=x.r_regionaware; R.loc[x.b,x.a]=x.r_regionaware
np.fill_diagonal(R.values,0.999); R.to_csv("coexpr_meanR_REGIONAWARE.csv")

# ---- impact on the FIGURE edges (POS r>=0.5 ; ANTI r<=-0.3 & padj<0.05) ----
POS,NEG=0.50,-0.30
def edges(rcol,pcol):
    pos=df[df[rcol]>=POS]; anti=df[(df[rcol]<=NEG)&(df[pcol]<0.05)]
    return set(map(tuple,pos[["a","b"]].values.tolist())), set(map(tuple,anti[["a","b"]].values.tolist()))
pos_ra,anti_ra=edges("r_regionaware","padj_regionaware")
pos_pl,anti_pl=edges("r_pooled","padj_pooled")
print("POS edges: region-aware=%d  pooled=%d"%(len(pos_ra),len(pos_pl)))
print("ANTI edges: region-aware=%d  pooled=%d"%(len(anti_ra),len(anti_pl)))
print("\nPOS lost when region-aware (pooled-only, likely region-driven):",sorted(pos_pl-pos_ra))
print("\nPOS gained by region-aware:",sorted(pos_ra-pos_pl))
print("\nANTI lost:",sorted(anti_pl-anti_ra)," ANTI gained:",sorted(anti_ra-anti_pl))
# DAB2 specifically
print("\n=== DAB2 pairs ===")
d2=df[(df.a=="DAB2")|(df.b=="DAB2")].copy()
d2["other"]=np.where(d2.a=="DAB2",d2.b,d2.a)
for _,x in d2.sort_values("r_regionaware").iterrows():
    if x.r_regionaware<=-0.3 or x.r_pooled<=-0.3:
        print(f"  DAB2-{x.other:8s} r_RA={x.r_regionaware:+.2f} padj_RA={x.padj_regionaware:.3f} | r_pool={x.r_pooled:+.2f} padj_pool={x.padj_pooled:.3f}")
