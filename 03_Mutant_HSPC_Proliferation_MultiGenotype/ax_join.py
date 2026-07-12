import os,sys,warnings,numpy as np,pandas as pd
warnings.filterwarnings("ignore"); from scipy import stats
from sklearn.decomposition import PCA; from sklearn.cluster import MiniBatchKMeans
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
TMP="/tmp/asxl1"; RES=sys.argv[1]
panel=pd.read_csv(os.path.join(TMP,"panel.csv"),header=None)[0].tolist(); gi={g:i for i,g in enumerate(panel)}
Xs=[]; obs=[]
for geno,tp in [("WT","wk04"),("KO","wk04"),("WT","wk36"),("KO","wk36")]:
    Xs.append(np.load(os.path.join(TMP,f"{geno}_{tp}_X.npy"))); obs.append(pd.read_csv(os.path.join(TMP,f"{geno}_{tp}_obs.csv")))
X=np.vstack(Xs); obs=pd.concat(obs,ignore_index=True)
# cell-cycle regression (clustering only)
C=np.column_stack([np.ones(len(obs)),obs["S_score"],obs["G2M_score"]])
Xr=X-C@np.linalg.lstsq(C,X,rcond=None)[0]
Z=np.clip((Xr-Xr.mean(0))/(Xr.std(0)+1e-9),-10,10)
pcs=PCA(20,random_state=0,svd_solver="randomized").fit_transform(Z)
obs["cluster"]=MiniBatchKMeans(12,random_state=0,n_init=3,batch_size=2048).fit_predict(pcs).astype(str)
sets={"HSC":["Procr","Hlf","Mecom","Fgd5","Mllt3","Hoxa9","Meis1","Mpl","Ctnnal1"],
 "MkP":["Pf4","Vwf","Itga2b","Gp1bb","Plek"],"Ery":["Klf1","Gata1","Car1","Car2","Gypa"],
 "GMP":["Elane","Mpo","Prtn3","Ctsg","Ms4a3"],"Lymph":["Dntt","Il7r","Vpreb1","Cd79a"]}
sig={k:X[:,[gi[g] for g in gs if g in gi]].mean(1) for k,gs in sets.items()}
sd=pd.DataFrame(sig); sd["cluster"]=obs["cluster"].values; clm=sd.groupby("cluster").mean()
clz=(clm-clm.mean())/clm.std(); lab=clz.idxmax(1).to_dict()
obs["celltype"]=obs["cluster"].map(lab); obs["cycling"]=obs["phase"]!="G1"
prim=[c for c in clz.index if clz.loc[c].idxmax()=="HSC"]
obs["compartment"]=np.where(obs["cluster"].isin(prim),"LT-HSC+MPP1","Committed")
obs.to_csv(os.path.join(RES,"asxl1_10x_joint_obs.csv"),index=False)
P=obs[obs.compartment=="LT-HSC+MPP1"]
print("cluster->lineage:",lab,"\nprimitive clusters:",prim)
print("\nLT-HSC+MPP1 cells:\n",P.groupby(["timepoint","genotype"]).size().to_string())
def comp(d,lab):
    w=d[d.genotype=="WT"]; k=d[d.genotype=="KO"]
    chi2,p,_,_=stats.chi2_contingency([[w.cycling.sum(),(~w.cycling).sum()],[k.cycling.sum(),(~k.cycling).sum()]])
    u,pu=stats.mannwhitneyu(k.prolif,w.prolif)
    print(f"  [{lab}]  cycling WT={w.cycling.mean()*100:.1f}% KO={k.cycling.mean()*100:.1f}%  (Δ{ (k.cycling.mean()-w.cycling.mean())*100:+.1f})  chi2p={p:.2g} | prolif WT={w.prolif.mean():.3f} KO={k.prolif.mean():.3f} MWUp={pu:.2g}")
    return (lab,w.cycling.mean()*100,k.cycling.mean()*100,p,w.prolif.mean(),k.prolif.mean(),pu)
print("\nWithin LT-HSC+MPP1 (batch-free, WT+KO co-captured per timepoint):")
rows=[comp(P[P.timepoint=="wk04"],"wk04"),comp(P[P.timepoint=="wk36"],"wk36"),comp(P,"pooled")]
print("\nWhole LSK (all compartments):")
comp(obs[obs.timepoint=="wk04"],"wk04-all"); comp(obs[obs.timepoint=="wk36"],"wk36-all")
# replicate-aware: pseudobulk per (genotype x timepoint), paired across timepoints
pb=P.groupby(["genotype","timepoint"]).agg(cyc=("cycling","mean"),prolif=("prolif","mean")).reset_index()
print("\nPseudobulk (LT-HSC+MPP1) per genotype x timepoint:"); print(pb.round(3).to_string(index=False))
piv=pb.pivot(index="timepoint",columns="genotype",values="cyc")
t,pp=stats.ttest_rel(piv["KO"],piv["WT"]); print(f"\nPaired KO-vs-WT cycling across timepoints: mean Δ={ (piv['KO']-piv['WT']).mean()*100:+.1f}pp  paired-p={pp:.3g} (n=2 pairs)")
pd.DataFrame(rows,columns=["strata","cyc_WT","cyc_KO","chi2_p","prolif_WT","prolif_KO","MWU_p"]).to_csv(os.path.join(RES,"asxl1_10x_comparison.csv"),index=False)
pb.to_csv(os.path.join(RES,"asxl1_10x_pseudobulk.csv"),index=False)
# composition
comp_tab=pd.crosstab(obs["compartment"],[obs["genotype"],obs["timepoint"]],normalize="columns")*100
print("\n% in each compartment:\n",comp_tab.round(1).to_string())
print("\nDONE")
