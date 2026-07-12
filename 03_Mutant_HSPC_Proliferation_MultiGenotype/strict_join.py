import os,sys,glob,numpy as np,pandas as pd; from scipy import stats
from sklearn.decomposition import PCA; from sklearn.cluster import MiniBatchKMeans
TMP=sys.argv[1]; unit=sys.argv[2]  # unit col: biorep or mouse
panel=pd.read_csv(f"{TMP}/panel.csv",header=None)[0].tolist(); gi={g:i for i,g in enumerate(panel)}
libs=sorted(set(os.path.basename(f)[:-6] for f in glob.glob(f"{TMP}/*_X.npy")))
X=np.vstack([np.load(f"{TMP}/{l}_X.npy") for l in libs]).astype("float32")
obs=pd.concat([pd.read_csv(f"{TMP}/{l}_obs.csv") for l in libs],ignore_index=True)
Xo=X.copy()
# cell-cycle regression (clustering only) + finer clustering
Cc=np.column_stack([np.ones(len(obs)),obs["S_score"],obs["G2M_score"]]).astype("float32")
Xr=X-Cc@np.linalg.solve(Cc.T@Cc,Cc.T@X); Xr-=Xr.mean(0); Xr/=(Xr.std(0)+1e-9); np.clip(Xr,-10,10,out=Xr)
idx=np.random.default_rng(0).choice(len(Xr),min(15000,len(Xr)),replace=False)
pcs=PCA(20,random_state=0,svd_solver="randomized").fit(Xr[idx]).transform(Xr)
obs["cluster"]=MiniBatchKMeans(20,random_state=0,n_init=3,batch_size=4096).fit_predict(pcs).astype(str)
# per-cell signatures on ORIGINAL expression, z-scored across cells
def zmean(genes):
    cols=[gi[g] for g in genes if g in gi]; sub=Xo[:,cols]
    z=(sub-sub.mean(0))/(sub.std(0)+1e-9); return z.mean(1)
sig={"HSC":zmean(["Procr","Fgd5","Hlf","Mecom","Ctnnal1","Mllt3","Hoxa9","Meis1","Mpl"]),
     "Mk":zmean(["Pf4","Vwf","Itga2b","Gp1bb","Plek"]),"Ery":zmean(["Klf1","Gata1","Car1","Car2","Gypa"]),
     "GMP":zmean(["Elane","Mpo","Prtn3","Ctsg","Ms4a3"]),"Ly":zmean(["Dntt","Il7r","Vpreb1","Cd79a"]),
     "Flt3":zmean(["Flt3"])}
S=pd.DataFrame(sig); S["cluster"]=obs["cluster"].values
cm=S.groupby("cluster").mean(); cz=(cm-cm.mean())/cm.std()
obs["cycling"]=obs["phase"]!="G1"
cyc=obs.groupby("cluster")["cycling"].mean()*100; frac=obs["cluster"].value_counts(normalize=True)*100
# LT-HSC+MPP1 = HSC-high, all lineages low, Flt3 low
linmax=cz[["Mk","Ery","GMP","Ly"]].max(1)
prim=[c for c in cz.index if cz.loc[c,"HSC"]>0.5 and linmax[c]<0.5]
tab=pd.DataFrame({"HSC_z":cz["HSC"].round(2),"maxLin_z":linmax.round(2),"Flt3_z":cz["Flt3"].round(2),
                  "pct_cells":frac.round(1),"cycling%":cyc.round(0)}).sort_values("HSC_z",ascending=False)
print("=== clusters (sorted by HSC_z) ==="); print(tab.to_string())
print("\nLT-HSC+MPP1 clusters:",prim)
obs["primitive"]=obs["cluster"].isin(prim)
P=obs[obs.primitive]
print(f"LT-HSC+MPP1 = {obs.primitive.mean()*100:.1f}% of cells | cycling {P.cycling.mean()*100:.1f}%")
obs.to_csv(f"{TMP}/strict_obs.csv",index=False)
# mouse-level
pb=P.groupby([unit,"genotype"]+(["timepoint"] if "timepoint" in P else [])).agg(cyc=("cycling","mean"),prolif=("prolif","mean"),n=("cycling","size")).reset_index()
pb=pb[pb.n>=30]
print("\nper-unit (LT-HSC+MPP1):"); print(pb.round(3).to_string(index=False))
def tt(s,l):
    w=s[s.genotype.isin(["WT"])]; k=s[s.genotype.isin(["KO","Dnmt3a","Jak2V617F"])]
    if len(w)>1 and len(k)>1:
        _,pc=stats.ttest_ind(k.cyc,w.cyc); _,pp=stats.ttest_ind(k.prolif,w.prolif)
        print(f"[{l}] WT={w.cyc.mean()*100:.1f}% MUT={k.cyc.mean()*100:.1f}% p_cyc={pc:.3g} | prolif p={pp:.3g} (n {len(w)}v{len(k)})")
if "timepoint" in pb:
    for tp in sorted(pb.timepoint.unique()): tt(pb[pb.timepoint==tp],tp)
tt(pb,"pooled")
pb.to_csv(f"{TMP}/strict_pb.csv",index=False)
