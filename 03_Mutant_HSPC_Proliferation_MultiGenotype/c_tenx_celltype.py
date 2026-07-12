import os,sys,gzip,tarfile,tempfile,shutil,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); RES=sys.argv[1]; geno=sys.argv[2]; pref=sys.argv[3]
td=tempfile.mkdtemp()
with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
    for suf in ["matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz"]:
        open(os.path.join(td,suf[:-3]),"wb").write(gzip.open(tf.extractfile(pref+suf)).read())
mtx=sio.mmread(os.path.join(td,"matrix.mtx")).tocsr()
feats=[l.split("\t")[1] for l in open(os.path.join(td,"features.tsv")).read().splitlines()]
bcs=open(os.path.join(td,"barcodes.tsv")).read().splitlines(); shutil.rmtree(td)
ad=sc.AnnData(mtx.T.tocsr().astype("float32")); ad.var_names=feats; ad.obs_names=bcs; ad.var_names_make_unique()
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None)
ad=ad[ad.obs["pct_counts_mt"]<10].copy(); sc.pp.filter_cells(ad,min_genes=300)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif_core")
markers={"HSC":["Procr","Hlf","Mecom","Fgd5","Slamf1","Mllt3","Hoxa9","Gata2"],
 "MPP":["Flt3","Cd34","Cd48"],"MkP":["Pf4","Vwf","Itga2b","Gp1bb","Plek"],
 "Ery":["Klf1","Gata1","Car1","Car2","Hba-a1","Gypa"],
 "GMP_Mye":["Elane","Mpo","Prtn3","Ctsg","Ms4a3","Cebpe"],
 "Lymph":["Dntt","Il7r","Vpreb1","Ebf1","Cd79a"],"Baso_MC":["Prss34","Mcpt8","Cma1"]}
for k,gs in markers.items(): sc.tl.score_genes(ad,[g for g in gs if g in ad.var_names],score_name="sc_"+k)
S=ad.obs[["sc_"+k for k in markers]]; ad.obs["celltype"]=S.idxmax(axis=1).str.replace("sc_","")
ad.obs["genotype"]=geno; ad.obs["cycling"]=(ad.obs["phase"]!="G1")
ad.obs[["genotype","celltype","phase","cycling","prolif_core"]].to_csv(os.path.join(RES,f"tenx_celltype_{geno}.csv"))
print(geno,"n=",ad.n_obs); print(ad.obs["celltype"].value_counts().to_string())
