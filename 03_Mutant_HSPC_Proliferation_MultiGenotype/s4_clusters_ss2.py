import os,sys,warnings,numpy as np,pandas as pd,scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,os.path.dirname(__file__))
from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(DATA,"results")
plates={"WT":["GSE124822_LTHSC-WT_476_RNA_counts.txt.gz","GSE124822_LTHSC-WT_846_RNA_counts.txt.gz"],
        "Tet2":["GSE124822_LTHSC-TET2_30_RNA_counts.txt.gz","GSE124822_LTHSC_TET2_272_RNA_counts.txt.gz"],
        "Dnmt3a":["GSE124822_LTHSC-DNMT3A_675_RNA_counts.txt.gz"]}
ads=[]
for geno,files in plates.items():
    mats=[pd.read_csv(os.path.join(DATA,f),sep="\t",index_col=0) for f in files]
    mats=[m.groupby(m.index).sum() for m in mats]
    genes=sorted(set().union(*[set(m.index) for m in mats]))
    expr=pd.concat([m.reindex(genes).fillna(0) for m in mats],axis=1)
    a=sc.AnnData(expr.T.astype("float32")); a.var_names_make_unique(); a.obs["genotype"]=geno
    a.obs_names=[f"{geno}_{i}_{b}" for i,b in enumerate(a.obs_names)]
    ads.append(a)
ad=sc.concat(ads,join="inner"); sc.pp.filter_cells(ad,min_genes=500)
ad.layers["counts"]=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); ad.raw=ad
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif_core")
sc.pp.highly_variable_genes(ad,n_top_genes=2000); adh=ad[:,ad.var.highly_variable].copy()
sc.pp.scale(adh,max_value=10); sc.tl.pca(adh,n_comps=30)
sc.pp.neighbors(adh,n_neighbors=15); sc.tl.leiden(adh,resolution=0.6); sc.tl.umap(adh)
ad.obs["leiden"]=adh.obs["leiden"]; ad.obsm["X_umap"]=adh.obsm["X_umap"]
# marker lineage scores
markers={"HSC":["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Ly6a"],
         "MPP":["Flt3","Cd34","Cd48"],"MkP":["Pf4","Vwf","Itga2b","Gp1bb"],
         "Ery":["Klf1","Gata1","Car1","Car2","Hba-a1"],"Mye":["Elane","Mpo","Prtn3","Ctsg","Gfi1"],
         "Lymph":["Dntt","Il7r","Vpreb1","Ebf1"]}
for k,gs in markers.items():
    sc.tl.score_genes(ad,[g for g in gs if g in ad.var_names],score_name="sc_"+k)
clm=ad.obs.groupby("leiden")[["sc_"+k for k in markers]].mean()
lab={cl:clm.loc[cl].idxmax().replace("sc_","") for cl in clm.index}
ad.obs["celltype"]=ad.obs["leiden"].map(lab)
# per cluster x genotype
ad.obs["cycling"]=(ad.obs["phase"]!="G1")
tab=ad.obs.groupby(["celltype","genotype"]).agg(n=("cycling","size"),
      cycling_pct=("cycling",lambda x:100*x.mean()),prolif=("prolif_core","mean")).reset_index()
tab.to_csv(os.path.join(RES,"ss2_cluster_proliferation.csv"),index=False)
clm.to_csv(os.path.join(RES,"ss2_cluster_marker_scores.csv"))
print("cluster labels:",lab)
print(tab.to_string(index=False))
print("\ncells per cluster:\n",ad.obs["celltype"].value_counts().to_string())
# figures
fig,axs=plt.subplots(1,3,figsize=(15,4.3))
sc.pl.umap(ad,color="genotype",ax=axs[0],show=False,title="genotype",frameon=False)
sc.pl.umap(ad,color="celltype",ax=axs[1],show=False,title="cluster (marker-annotated)",frameon=False)
sc.pl.umap(ad,color="phase",ax=axs[2],show=False,title="cell-cycle phase",frameon=False)
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig4_ss2_umap.png"),bbox_inches="tight",dpi=130)
ad.write(os.path.join(RES,"ss2_clustered.h5ad"))
print("done")
