import os,sys,gzip,tarfile,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); TMP="/tmp/dnmt3a"
gsm,geno,sex,rep,mode=sys.argv[1:6]
MARK=["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Gata2","Meis1","Mpl","Ctnnal1",
 "Flt3","Cd34","Cd48","Mki67","Pf4","Vwf","Itga2b","Gp1bb","Plek","Klf1","Gata1","Car1","Car2","Gypa","Epor",
 "Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe","Dntt","Il7r","Vpreb1","Ebf1","Cd79a","Prss34","Mcpt8","Cpa3"]
with tarfile.open(os.path.join(DATA,"GSE272266_RAW.tar")) as tf:
    mem={m.name.split("_")[-2]+"_"+m.name.split("_")[-1].split(".")[0]:m for m in tf.getmembers() if m.name.startswith(gsm)}
    def g(k): return [m for n,m in mem.items() if k in n][0]
    M=sio.mmread(gzip.open(tf.extractfile(g("matrix")))).tocsr()
    feats=[l.split("\t")[1] for l in gzip.open(tf.extractfile(g("features")),"rt").read().splitlines()]
    bcs=gzip.open(tf.extractfile(g("barcodes")),"rt").read().splitlines()
ad=sc.AnnData(M.T.tocsr().astype("float32")); ad.var_names=feats; ad.var_names_make_unique()
ad.obs_names=[f"{geno}_{sex}{rep}_{b}" for b in bcs]
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[(ad.obs["pct_counts_mt"]<10)].copy(); sc.pp.filter_cells(ad,min_genes=500)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[x for x in S_PHASE if x in ad.var_names]; g2=[x for x in G2M_PHASE if x in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[x for x in PROLIFERATION_CORE if x in ad.var_names],score_name="prolif")
if mode=="define":
    sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor="seurat")
    panel=sorted(set(ad.var_names[ad.var.highly_variable]).union([m for m in MARK if m in ad.var_names]))
    pd.Series(panel).to_csv(os.path.join(TMP,"panel.csv"),index=False,header=False)
else:
    panel=pd.read_csv(os.path.join(TMP,"panel.csv"),header=None)[0].tolist()
np.save(os.path.join(TMP,f"{geno}_{sex}{rep}_X.npy"),np.asarray(ad[:,panel].X.todense(),dtype="float32"))
ad.obs[["S_score","G2M_score","phase","prolif"]].assign(genotype=geno,sex=sex,rep=rep,barcode=ad.obs_names).to_csv(os.path.join(TMP,f"{geno}_{sex}{rep}_obs.csv"),index=False)
print(gsm,geno,sex,rep,mode,"cells",ad.n_obs)
