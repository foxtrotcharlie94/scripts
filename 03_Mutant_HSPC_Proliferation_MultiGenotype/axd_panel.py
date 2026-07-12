import os,sys,glob,gzip,re,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); TMP="/tmp/axd"
tag,mode=sys.argv[1],sys.argv[2]   # tag e.g. wk04_wt
MARK=["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Gata2","Meis1","Mpl","Ctnnal1","Flt3","Cd34","Cd48","Mki67",
 "Pf4","Vwf","Itga2b","Gp1bb","Plek","Klf1","Gata1","Car1","Car2","Gypa","Epor","Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe",
 "Dntt","Il7r","Vpreb1","Ebf1","Cd79a","Prss34","Mcpt8","Cpa3"]
def f(k): return glob.glob(os.path.join(DATA,f"*{tag}*{k}*"))[0]
M=sio.mmread(gzip.open(f("matrix"))).tocsr()
fl=[l.split("\t") for l in gzip.open(f("features"),"rt").read().splitlines()]
ftype=np.array([c[2] if len(c)>2 else "Gene Expression" for c in fl]); fsym=np.array([c[1] for c in fl])
bcs=gzip.open(f("barcodes"),"rt").read().splitlines()
gex=np.where(ftype=="Gene Expression")[0]
hidx=[i for i,n in enumerate(fsym) if re.match(r"^wk\d+_(wt|mut)_\d+$",n)]   # sample hashtags only
hnames=fsym[hidx]
H=np.asarray(M[hidx,:].todense())                    # nHTO x cells
top=H.argmax(0); tot=H.sum(0); frac=np.where(tot>0,H.max(0)/(tot+1e-9),0)
biorep=np.array([hnames[t] for t in top]); geno=np.where(np.char.find(biorep.astype(str),"_wt_")>=0,"WT","KO")
ad=sc.AnnData(M[gex,:].T.tocsr().astype("float32")); ad.var_names=list(fsym[gex]); ad.var_names_make_unique()
ad.obs["biorep"]=biorep; ad.obs["genotype"]=geno; ad.obs["hto_frac"]=frac; ad.obs["hto_tot"]=tot
ad.obs_names=[f"{b}_{bc}" for b,bc in zip(biorep,bcs)]
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[(ad.obs.pct_counts_mt<10)&(ad.obs.hto_frac>0.6)&(ad.obs.hto_tot>=3)].copy(); sc.pp.filter_cells(ad,min_genes=300)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[x for x in S_PHASE if x in ad.var_names]; g2=[x for x in G2M_PHASE if x in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2); sc.tl.score_genes(ad,[x for x in PROLIFERATION_CORE if x in ad.var_names],score_name="prolif")
tp=tag.split("_")[0]
if mode=="define":
    sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor="seurat")
    panel=sorted(set(ad.var_names[ad.var.highly_variable]).union([m for m in MARK if m in ad.var_names]))
    pd.Series(panel).to_csv(os.path.join(TMP,"panel.csv"),index=False,header=False)
else:
    panel=pd.read_csv(os.path.join(TMP,"panel.csv"),header=None)[0].tolist()
np.save(os.path.join(TMP,f"{tag}_X.npy"),np.asarray(ad[:,panel].X.todense(),dtype="float32"))
ad.obs[["biorep","genotype","phase","prolif","S_score","G2M_score"]].assign(timepoint=tp).to_csv(os.path.join(TMP,f"{tag}_obs.csv"),index=False)
print(tag,mode,"cells",ad.n_obs,"\n",ad.obs["biorep"].value_counts().to_string())
