import os,sys,warnings,numpy as np,pandas as pd,scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); RES=sys.argv[1]
plates={"WT":["GSE124822_LTHSC-WT_476_RNA_counts.txt.gz","GSE124822_LTHSC-WT_846_RNA_counts.txt.gz"],
        "Tet2":["GSE124822_LTHSC-TET2_30_RNA_counts.txt.gz","GSE124822_LTHSC_TET2_272_RNA_counts.txt.gz"],
        "Dnmt3a":["GSE124822_LTHSC-DNMT3A_675_RNA_counts.txt.gz"]}
ads=[]
for geno,files in plates.items():
    mats=[pd.read_csv(os.path.join(DATA,f),sep="\t",index_col=0) for f in files]
    mats=[m.groupby(m.index).sum() for m in mats]
    g=sorted(set().union(*[set(m.index) for m in mats]))
    expr=pd.concat([m.reindex(g).fillna(0) for m in mats],axis=1)
    a=sc.AnnData(expr.T.astype("float32")); a.var_names_make_unique(); a.obs["genotype"]=geno
    a.obs_names=[f"{geno}_{i}" for i in range(a.n_obs)]; ads.append(a)
ad=sc.concat(ads,join="inner"); sc.pp.filter_cells(ad,min_genes=500)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[x for x in S_PHASE if x in ad.var_names]; g2=[x for x in G2M_PHASE if x in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[x for x in PROLIFERATION_CORE if x in ad.var_names],score_name="prolif_core")
markers={"HSC":["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Ly6a"],"MPP":["Flt3","Cd34"],
 "MkP":["Pf4","Vwf","Itga2b","Gp1bb"],"Ery":["Klf1","Gata1","Car1","Car2","Hba-a1"],
 "Mye":["Elane","Mpo","Prtn3","Ctsg","Gfi1"],"Lymph":["Dntt","Il7r","Vpreb1","Ebf1"]}
for k,gs in markers.items(): sc.tl.score_genes(ad,[x for x in gs if x in ad.var_names],score_name="sc_"+k)
sc.pp.highly_variable_genes(ad,n_top_genes=2000); adh=ad[:,ad.var.highly_variable].copy()
sc.pp.scale(adh,max_value=10); sc.tl.pca(adh,n_comps=30); ad.obsm["X_pca"]=adh.obsm["X_pca"]
# KMeans clustering on PCA (fast, no numba)
from sklearn.cluster import KMeans
ad.obs["leiden"]=KMeans(n_clusters=6,random_state=0,n_init=10).fit_predict(ad.obsm["X_pca"]).astype(str)
mk=list(markers); clm=ad.obs.groupby("leiden")[["sc_"+k for k in mk]].mean()
seen={}; lab={}
for cl in clm.index:
    b=clm.loc[cl].idxmax().replace("sc_",""); seen[b]=seen.get(b,0)+1; lab[cl]=b if seen[b]==1 else f"{b}{seen[b]}"
ad.obs["celltype"]=ad.obs["leiden"].map(lab); ad.obs["cycling"]=(ad.obs["phase"]!="G1")
tab=ad.obs.groupby(["celltype","genotype"]).agg(n=("cycling","size"),
     cycling_pct=("cycling",lambda x:round(100*x.mean(),1)),prolif=("prolif_core",lambda x:round(x.mean(),3))).reset_index()
tab.to_csv(os.path.join(RES,"ss2_cluster_proliferation.csv"),index=False)
clm.round(3).to_csv(os.path.join(RES,"ss2_cluster_marker_scores.csv"))
ad.obs[["genotype","leiden","celltype","phase","prolif_core"]].to_csv(os.path.join(RES,"ss2_obs_clustered.csv"))
print("cluster->lineage:",lab); print()
print("cells/cluster x genotype:"); print(pd.crosstab(ad.obs.celltype,ad.obs.genotype).to_string()); print()
print(tab.to_string(index=False))
# PCA scatter figure
pca=ad.obsm["X_pca"]
fig,axs=plt.subplots(1,3,figsize=(15,4.3))
for ax,key,title in zip(axs,["genotype","celltype","phase"],["genotype","cluster","phase"]):
    cats=ad.obs[key].astype("category"); 
    for c in cats.cat.categories:
        m=(cats==c).values; ax.scatter(pca[m,0],pca[m,1],s=14,label=str(c),alpha=0.8)
    ax.set_title(title); ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(fontsize=7,markerscale=1.5)
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig4_ss2_pca.png"),bbox_inches="tight",dpi=130)
print("FIG_OK")
