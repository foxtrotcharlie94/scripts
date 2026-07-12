import os,sys,warnings,numpy as np,pandas as pd
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA; from sklearn.cluster import KMeans
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
TMP="/tmp/hspc"; RES=sys.argv[1]
panel=pd.read_csv(os.path.join(TMP,"tenx_panel.csv"),header=None)[0].tolist()
Xw=np.load(os.path.join(TMP,"tenx_WT_panel.npy")); Xv=np.load(os.path.join(TMP,"tenx_Jak2V617F_panel.npy"))
ow=pd.read_csv(os.path.join(TMP,"tenx_WT_obs.csv")); ov=pd.read_csv(os.path.join(TMP,"tenx_Jak2V617F_obs.csv"))
X=np.vstack([Xw,Xv]); obs=pd.concat([ow,ov],ignore_index=True)
gi={g:i for i,g in enumerate(panel)}
# ---- cell-cycle regression (clustering only) ----
C=np.column_stack([np.ones(len(obs)),obs["S_score"].values,obs["G2M_score"].values])
beta,_,_,_=np.linalg.lstsq(C,X,rcond=None); Xr=X-C@beta
# scale (z per gene) then PCA
mu=Xr.mean(0); sd=Xr.std(0); sd[sd==0]=1; Xs=np.clip((Xr-mu)/sd,-10,10)
pcs=PCA(n_components=30,random_state=0).fit_transform(Xs)
obs["cluster"]=KMeans(n_clusters=12,random_state=0,n_init=10).fit_predict(pcs).astype(str)
# ---- annotate clusters by lineage signatures (original log-norm panel X) ----
sets={"HSC":["Procr","Hlf","Mecom","Fgd5","Mllt3","Hoxa9","Meis1","Mpl","Ctnnal1"],
 "MkP":["Pf4","Vwf","Itga2b","Gp1bb","Plek"],"Ery":["Klf1","Gata1","Car1","Car2","Gypa","Epor"],
 "GMP":["Elane","Mpo","Prtn3","Ctsg","Ms4a3"],"Lymph":["Dntt","Il7r","Vpreb1","Cd79a"],
 "Baso":["Prss34","Mcpt8","Cpa3"]}
sig={}
for k,gs in sets.items():
    cols=[gi[g] for g in gs if g in gi]; sig[k]=X[:,cols].mean(1)
sigdf=pd.DataFrame(sig); sigdf["cluster"]=obs["cluster"].values
clmean=sigdf.groupby("cluster").mean()
clz=(clmean-clmean.mean())/clmean.std()   # z across clusters
lab=clz.idxmax(axis=1).to_dict()
obs["celltype"]=obs["cluster"].map(lab)
obs["Flt3"]=X[:,gi["Flt3"]] if "Flt3" in gi else 0
obs["cycling"]=obs["phase"]!="G1"
# LT-HSC+MPP1 = clusters whose dominant signature is HSC
prim=[c for c in clz.index if clz.loc[c].idxmax()=="HSC"]
obs["compartment"]=np.where(obs["cluster"].isin(prim),"LT-HSC+MPP1","Committed/other")
obs.to_csv(os.path.join(RES,"tenx_joint_obs.csv"),index=False)
clz.round(2).to_csv(os.path.join(RES,"tenx_cluster_signatures.csv"))
# ---- composition + proliferation within compartment ----
comp=pd.crosstab(obs["compartment"],obs["genotype"],normalize="columns")*100
print("cluster->lineage:",lab); print("primitive clusters:",prim)
print("\n% of each genotype in compartment:\n",comp.round(1).to_string())
def cmp(df,label):
    w=df[df.genotype=="WT"]; v=df[df.genotype=="Jak2V617F"]
    tab=np.array([[w.cycling.sum(),(~w.cycling).sum()],[v.cycling.sum(),(~v.cycling).sum()]])
    chi2,p,_,_=stats.chi2_contingency(tab); u,pu=stats.mannwhitneyu(v.prolif_core,w.prolif_core)
    print(f"\n[{label}] n_WT={len(w)} n_VF={len(v)}  cycling WT={w.cycling.mean()*100:.1f}% VF={v.cycling.mean()*100:.1f}%  chi2_p={p:.2e}"
          f"  prolif WT={w.prolif_core.mean():.3f} VF={v.prolif_core.mean():.3f} MWU_p={pu:.2e}")
    return dict(compartment=label,n_WT=len(w),n_VF=len(v),cycling_WT=round(w.cycling.mean()*100,1),
        cycling_VF=round(v.cycling.mean()*100,1),chi2_p=p,prolif_WT=round(w.prolif_core.mean(),3),
        prolif_VF=round(v.prolif_core.mean(),3),MWU_p=pu)
res=[cmp(obs[obs.compartment=="LT-HSC+MPP1"],"LT-HSC+MPP1"),
     cmp(obs[obs.compartment=="Committed/other"],"Committed/other"),
     cmp(obs,"All LSK")]
pd.DataFrame(res).to_csv(os.path.join(RES,"tenx_LTHSC_MPP1_comparison.csv"),index=False)
# ---- figure ----
fig,axs=plt.subplots(1,3,figsize=(16,4.6))
for c in sorted(set(obs["celltype"])):
    m=(obs["celltype"]==c).values; axs[0].scatter(pcs[m,0],pcs[m,1],s=4,alpha=.5,label=c)
axs[0].set_title("clusters (CC-regressed)"); axs[0].legend(fontsize=7,markerscale=2); axs[0].set_xlabel("PC1");axs[0].set_ylabel("PC2")
for c,col in [("WT","#7f8c8d"),("Jak2V617F","#c0392b")]:
    m=(obs["genotype"]==c).values; axs[1].scatter(pcs[m,0],pcs[m,1],s=4,alpha=.4,label=c,color=col)
axs[1].set_title("genotype"); axs[1].legend(markerscale=2)
prim_mask=obs["compartment"].eq("LT-HSC+MPP1").values
axs[2].scatter(pcs[~prim_mask,0],pcs[~prim_mask,1],s=4,alpha=.3,color="#d0d0d0",label="committed/other")
axs[2].scatter(pcs[prim_mask,0],pcs[prim_mask,1],s=4,alpha=.5,color="#2c7fb8",label="LT-HSC+MPP1")
axs[2].set_title("LT-HSC+MPP1 compartment"); axs[2].legend(markerscale=2)
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig5_tenx_clusters.png"),bbox_inches="tight",dpi=130)
print("\nFIG_OK")
