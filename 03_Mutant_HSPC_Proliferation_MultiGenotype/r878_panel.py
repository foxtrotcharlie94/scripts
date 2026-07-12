import os,sys,tarfile,tempfile,gzip,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); TMP="/tmp/r878"
gsm,geno,mouse,mode=sys.argv[1:5]
MARK=["Procr","Fgd5","Hlf","Mecom","Ctnnal1","Mllt3","Hoxa9","Meis1","Mpl","Slamf1","Gata2","Cd34","Flt3","Cd48","Mki67",
 "Klf1","Gata1","Car1","Car2","Gypa","Epor","Pf4","Vwf","Itga2b","Gp1bb","Plek","Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe",
 "Ltf","Camp","Ngp","S100a8","S100a9","Retnlg","Csf1r","Ly6c2","F13a1","Irf8","Siglech","Cd74","H2-Aa",
 "Prss34","Mcpt8","Cpa3","Cd200r3","Prg2","Dntt","Vpreb1","Vpreb3","Ebf1","Il7r","Rag1","Cd79a","Cd79b","Cd19","Ms4a1","Ighm",
 "Cd3e","Cd3d","Cd5","Cd6","Trac","Cd8a","Nkg7","Klrb1c","Ncr1","Gzmb","Klrd1"]
td=tempfile.mkdtemp()
with tarfile.open(os.path.join(DATA,f"{gsm}.tar")) as tf:
    for m in tf.getmembers(): tf.extract(m,td)
import glob
mtx=glob.glob(f"{td}/*matrix.mtx.gz")[0]; feat=glob.glob(f"{td}/*features.tsv.gz")[0]; bc=glob.glob(f"{td}/*barcodes.tsv.gz")[0]
M=sio.mmread(gzip.open(mtx)).tocsr()
fl=[l.split("\t") for l in gzip.open(feat,"rt").read().splitlines()]
ftype=np.array([c[2] if len(c)>2 else "Gene Expression" for c in fl]); fsym=np.array([c[1] for c in fl])
bcs=gzip.open(bc,"rt").read().splitlines()
keep=np.where(ftype=="Gene Expression")[0]; M=M[keep,:]; feats=list(fsym[keep])
ad=sc.AnnData(M.T.tocsr().astype("float32")); ad.var_names=feats; ad.var_names_make_unique()
ad.obs_names=[f"{mouse}_{b}" for b in bcs]; ad.obs["mouse"]=mouse; ad.obs["genotype"]=geno
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[(ad.obs.pct_counts_mt<10)].copy(); sc.pp.filter_cells(ad,min_genes=500)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2); sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif")
if mode=="define":
    sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor="seurat")
    panel=sorted(set(ad.var_names[ad.var.highly_variable]).union([m for m in MARK if m in ad.var_names]))
    pd.Series(panel).to_csv(f"{TMP}/panel.csv",index=False,header=False)
else: panel=pd.read_csv(f"{TMP}/panel.csv",header=None)[0].tolist()
np.save(f"{TMP}/{mouse}_X.npy",np.asarray(ad[:,panel].X.todense(),dtype="float32"))
ad.obs[["mouse","genotype","phase","prolif","S_score","G2M_score"]].to_csv(f"{TMP}/{mouse}_obs.csv",index=False)
print(gsm,geno,mouse,mode,"cells",ad.n_obs)
