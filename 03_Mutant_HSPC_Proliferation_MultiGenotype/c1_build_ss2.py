import os,sys,warnings,numpy as np,pandas as pd,scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
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
    a.obs_names=[f"{geno}_{i}" for i in range(a.n_obs)]; ads.append(a)
ad=sc.concat(ads,join="inner"); sc.pp.filter_cells(ad,min_genes=500)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif_core")
markers={"HSC":["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Ly6a"],
         "MPP":["Flt3","Cd34"],"MkP":["Pf4","Vwf","Itga2b","Gp1bb"],
         "Ery":["Klf1","Gata1","Car1","Car2","Hba-a1"],"Mye":["Elane","Mpo","Prtn3","Ctsg","Gfi1"],
         "Lymph":["Dntt","Il7r","Vpreb1","Ebf1"]}
for k,gs in markers.items(): sc.tl.score_genes(ad,[g for g in gs if g in ad.var_names],score_name="sc_"+k)
sc.pp.highly_variable_genes(ad,n_top_genes=2000)
adh=ad[:,ad.var.highly_variable].copy(); sc.pp.scale(adh,max_value=10); sc.tl.pca(adh,n_comps=30)
ad.obsm["X_pca"]=adh.obsm["X_pca"]
ad.write(os.path.join(RES,"ss2_built.h5ad"))
print("built", ad.shape, "saved")
