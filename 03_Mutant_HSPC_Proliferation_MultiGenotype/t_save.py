import os,sys,gzip,tarfile,tempfile,shutil,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio, sparse
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); W=sys.argv[1]; geno=sys.argv[2]; pref=sys.argv[3]
td=tempfile.mkdtemp()
with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
    for suf in ["matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz"]:
        open(os.path.join(td,suf[:-3]),"wb").write(gzip.open(tf.extractfile(pref+suf)).read())
mtx=sio.mmread(os.path.join(td,"matrix.mtx")).tocsr()
feats=[l.split("\t")[1] for l in open(os.path.join(td,"features.tsv")).read().splitlines()]
bcs=open(os.path.join(td,"barcodes.tsv")).read().splitlines(); shutil.rmtree(td)
ad=sc.AnnData(mtx.T.tocsr().astype("float32")); ad.var_names=feats; ad.obs_names=[f"{geno}_{b}" for b in bcs]; ad.var_names_make_unique()
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[ad.obs["pct_counts_mt"]<10].copy(); sc.pp.filter_cells(ad,min_genes=300)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif_core")
sparse.save_npz(os.path.join(W,f"tenx_{geno}_lognorm.npz"), ad.X.tocsr())
ad.obs[["genotype","S_score","G2M_score","phase","prolif_core"]].assign(genotype=geno).to_csv(os.path.join(W,f"tenx_{geno}_obs.csv"))
pd.Series(ad.var_names).to_csv(os.path.join(W,"tenx_varnames.csv"),index=False,header=False)
print(geno,"saved",ad.shape)
