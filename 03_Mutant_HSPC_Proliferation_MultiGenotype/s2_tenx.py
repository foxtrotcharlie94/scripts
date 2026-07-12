import os,sys,gzip,tarfile,tempfile,shutil,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.path.dirname(__file__))
from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(DATA,"results")
geno,pref=sys.argv[1],sys.argv[2]
td=tempfile.mkdtemp()
with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
    for suf in ["matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz"]:
        open(os.path.join(td,suf[:-3]),"wb").write(gzip.open(tf.extractfile(pref+suf)).read())
print("decompressed",flush=True)
mtx=sio.mmread(os.path.join(td,"matrix.mtx")).tocsr()
feats=[l.split("\t")[1] for l in open(os.path.join(td,"features.tsv")).read().splitlines()]
bcs=open(os.path.join(td,"barcodes.tsv")).read().splitlines(); shutil.rmtree(td)
print("mmread",mtx.shape,flush=True)
ad=sc.AnnData(mtx.T.tocsr().astype("float32")); ad.var_names=feats; ad.obs_names=bcs; ad.var_names_make_unique()
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[ad.obs["pct_counts_mt"]<10].copy(); sc.pp.filter_cells(ad,min_genes=300)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
pc=[g for g in PROLIFERATION_CORE if g in ad.var_names]; sc.tl.score_genes(ad,pc,score_name="prolif_core")
ad.obs[["phase","S_score","G2M_score","prolif_core"]].to_csv(os.path.join(RES,f"obs_tenx_{geno}.csv"))
ph=ad.obs["phase"].value_counts(normalize=True)*100
print(geno,"n=",ad.n_obs,"G1%%=%.1f S%%=%.1f G2M%%=%.1f cycling%%=%.1f prolif=%.3f"%(
    ph.get("G1",0),ph.get("S",0),ph.get("G2M",0),(ad.obs["phase"]!="G1").mean()*100,ad.obs["prolif_core"].mean()),flush=True)
