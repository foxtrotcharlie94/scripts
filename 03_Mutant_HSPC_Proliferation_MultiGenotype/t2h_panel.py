import os,sys,glob,gzip,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); TMP="/tmp/tet2h"
group,lane,geno,cond,mode=sys.argv[1:6]   # e.g. WT_Unt rep1_rep2 WT Unt define
MARK=["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Gata2","Meis1","Mpl","Ctnnal1",
 "Flt3","Cd34","Cd48","Mki67","Pf4","Vwf","Itga2b","Gp1bb","Plek","Klf1","Gata1","Car1","Car2","Gypa","Epor",
 "Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe","Dntt","Il7r","Vpreb1","Ebf1","Cd79a","Prss34","Mcpt8","Cpa3"]
pre=f"GSE209994_{group}_{lane}_filtered_feature_bc_matrix"
def f(k): return glob.glob(os.path.join(DATA,f"{pre}_{k}*"))[0]
M=sio.mmread(gzip.open(f("matrix"))).tocsr()                      # features x cells
fl=[l.split("\t") for l in gzip.open(f("features"),"rt").read().splitlines()]
ftype=np.array([c[2] if len(c)>2 else "Gene Expression" for c in fl]); fsym=np.array([c[1] for c in fl])
bcs=gzip.open(f("barcodes"),"rt").read().splitlines()
gex=np.where(ftype=="Gene Expression")[0]; hto=np.where(ftype=="Antibody Capture")[0]
H=np.asarray(M[hto,:].todense())                                 # 2 x cells
assign=H.argmax(0); frac=H.max(0)/(H.sum(0)+1e-9)
ad=sc.AnnData(M[gex,:].T.tocsr().astype("float32")); ad.var_names=list(fsym[gex]); ad.var_names_make_unique()
# biorep: lane rep1_rep2 -> h0=rep1,h1=rep2 ; rep3_rep4 -> h0=rep3,h1=rep4
reps=lane.split("_"); repname=np.array([reps[a] for a in assign])
ad.obs["mouse"]=[f"{geno}_{cond}_{r}" for r in repname]; ad.obs["hto_frac"]=frac
ad.obs_names=[f"{geno}_{cond}_{r}_{b}" for r,b in zip(repname,bcs)]
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[(ad.obs.pct_counts_mt<10)&(ad.obs.hto_frac>0.6)].copy(); sc.pp.filter_cells(ad,min_genes=500)  # drop HTO doublets
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
np.save(os.path.join(TMP,f"{group}_{lane}_X.npy"),np.asarray(ad[:,panel].X.todense(),dtype="float32"))
ad.obs[["mouse","genotype" if False else "hto_frac","phase","prolif","S_score","G2M_score"]].assign(genotype=geno,cond=cond).to_csv(os.path.join(TMP,f"{group}_{lane}_obs.csv"),index=False)
print(group,lane,mode,"cells",ad.n_obs,"| mice:",dict(ad.obs["mouse"].value_counts()))
