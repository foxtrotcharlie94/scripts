import os,sys,tarfile,tempfile,warnings,numpy as np,pandas as pd,scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); TMP="/tmp/jak2"
gsm,geno,rep,mode=sys.argv[1:5]
MARK=["Procr","Fgd5","Hlf","Mecom","Ctnnal1","Mllt3","Hoxa9","Meis1","Mpl","Slamf1","Gata2","Flt3","Cd34","Cd48","Mki67",
 "Pf4","Vwf","Itga2b","Gp1bb","Plek","Klf1","Gata1","Car1","Car2","Gypa","Epor","Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe",
 "Dntt","Il7r","Vpreb1","Ebf1","Cd79a","Prss34","Mcpt8","Cpa3"]
td=tempfile.mkdtemp()
with tarfile.open(os.path.join(DATA,f"{gsm}.tar")) as tf:
    m=[x for x in tf.getmembers() if x.name.endswith(".h5")][0]; tf.extract(m,td); h5=os.path.join(td,m.name)
ad=sc.read_10x_h5(h5); ad.var_names_make_unique()
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[(ad.obs.pct_counts_mt<10)].copy(); sc.pp.filter_cells(ad,min_genes=1000)
ad.obs_names=[f"{geno}{rep}_{b}" for b in ad.obs_names]; ad.obs["mouse"]=f"Jak2{geno}_{rep}"; ad.obs["genotype"]="WT" if geno=="WT" else "MUT"
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2); sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif")
if mode=="define":
    sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor="seurat")
    panel=sorted(set(ad.var_names[ad.var.highly_variable]).union([m for m in MARK if m in ad.var_names]))
    pd.Series(panel).to_csv(f"{TMP}/panel.csv",index=False,header=False)
else: panel=pd.read_csv(f"{TMP}/panel.csv",header=None)[0].tolist()
np.save(f"{TMP}/{geno}{rep}_X.npy",np.asarray(ad[:,panel].X.todense(),dtype="float32"))
ad.obs[["mouse","genotype","phase","prolif","S_score","G2M_score"]].to_csv(f"{TMP}/{geno}{rep}_obs.csv",index=False)
print(gsm,geno,rep,mode,"cells",ad.n_obs)
