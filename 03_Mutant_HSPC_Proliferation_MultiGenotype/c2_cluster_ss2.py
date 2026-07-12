import os,warnings,numpy as np,pandas as pd,scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"results")
ad=sc.read_h5ad(os.path.join(RES,"ss2_built.h5ad"))
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=30)
sc.tl.leiden(ad,resolution=0.6,flavor="igraph",n_iterations=2,directed=False)
sc.tl.umap(ad)
mk=["HSC","MPP","MkP","Ery","Mye","Lymph"]
clm=ad.obs.groupby("leiden")[["sc_"+k for k in mk]].mean()
lab={cl:clm.loc[cl].idxmax().replace("sc_","") for cl in clm.index}
# disambiguate duplicate labels by index
seen={}; lab2={}
for cl in clm.index:
    base=lab[cl]; seen[base]=seen.get(base,0)+1
    lab2[cl]=base if seen[base]==1 else f"{base}{seen[base]}"
ad.obs["celltype"]=ad.obs["leiden"].map(lab2)
ad.obs["cycling"]=(ad.obs["phase"]!="G1")
tab=ad.obs.groupby(["celltype","genotype"]).agg(n=("cycling","size"),
     cycling_pct=("cycling",lambda x:round(100*x.mean(),1)),
     prolif=("prolif_core",lambda x:round(x.mean(),3))).reset_index()
tab.to_csv(os.path.join(RES,"ss2_cluster_proliferation.csv"),index=False)
clm.round(3).to_csv(os.path.join(RES,"ss2_cluster_marker_scores.csv"))
print("cluster -> lineage:",lab2)
print(ad.obs["celltype"].value_counts().to_string()); print()
print(tab.to_string(index=False))
fig,axs=plt.subplots(1,3,figsize=(15,4.3))
sc.pl.umap(ad,color="genotype",ax=axs[0],show=False,title="genotype",frameon=False,s=40)
sc.pl.umap(ad,color="celltype",ax=axs[1],show=False,title="cluster",frameon=False,s=40)
sc.pl.umap(ad,color="phase",ax=axs[2],show=False,title="phase",frameon=False,s=40)
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig4_ss2_umap.png"),bbox_inches="tight",dpi=130)
ad.write(os.path.join(RES,"ss2_clustered.h5ad")); print("done")
