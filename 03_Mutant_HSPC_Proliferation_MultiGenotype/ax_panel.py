import os,sys,glob,gzip,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); TMP="/tmp/asxl1"
tag,geno,tp,mode=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
MARK=["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Gata2","Meis1","Mpl","Ctnnal1",
 "Flt3","Cd34","Cd48","Mki67","Bmi1","Pf4","Vwf","Itga2b","Gp1bb","Plek",
 "Klf1","Gata1","Car1","Car2","Gypa","Epor","Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe",
 "Dntt","Il7r","Vpreb1","Ebf1","Cd79a","Prss34","Mcpt8","Cpa3"]
def find(kind): return glob.glob(os.path.join(DATA,f"*{tag}*{kind}*"))[0]
m=sio.mmread(gzip.open(find("matrix"))).tocsr()   # features x cells
flines=[l.split("\t") for l in gzip.open(find("features"),"rt").read().splitlines()]
ftype=np.array([(c[2] if len(c)>2 else "Gene Expression") for c in flines])
fsym=np.array([c[1] for c in flines])
keep=np.where(ftype=="Gene Expression")[0]
m=m[keep,:]; feats=list(fsym[keep])
bcs=gzip.open(find("barcodes"),"rt").read().splitlines()
ad=sc.AnnData(m.T.tocsr().astype("float32")); ad.var_names=feats; ad.var_names_make_unique()
ad.obs_names=[f"{geno}_{tp}_{b}" for b in bcs]
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[ad.obs["pct_counts_mt"]<10].copy(); sc.pp.filter_cells(ad,min_genes=300)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif")
if mode=="define":
    sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor="seurat")
    panel=sorted(set(ad.var_names[ad.var.highly_variable]).union([x for x in MARK if x in ad.var_names]))
    pd.Series(panel).to_csv(os.path.join(TMP,"panel.csv"),index=False,header=False)
else:
    panel=pd.read_csv(os.path.join(TMP,"panel.csv"),header=None)[0].tolist()
np.save(os.path.join(TMP,f"{geno}_{tp}_X.npy"),np.asarray(ad[:,panel].X.todense(),dtype="float32"))
ad.obs[["S_score","G2M_score","phase","prolif"]].assign(genotype=geno,timepoint=tp,barcode=ad.obs_names).to_csv(os.path.join(TMP,f"{geno}_{tp}_obs.csv"),index=False)
print(geno,tp,mode,"cells",ad.n_obs)
