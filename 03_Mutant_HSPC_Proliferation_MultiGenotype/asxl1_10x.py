# Batch-free single-cell Asxl1 WT-vs-KO proliferation (GSE158184 10x, HTO co-captured).
# Expects in the data folder (any subdir): for each of wk04_wt, wk04_mut, wk36_wt, wk36_mut
#   *<tag>*matrix.mtx.gz, *<tag>*barcodes.tsv.gz, *<tag>*features.tsv.gz
import os,sys,glob,gzip,warnings,numpy as np,pandas as pd,scanpy as sc
from scipy import io as sio, stats
from sklearn.decomposition import PCA; from sklearn.cluster import KMeans
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.path.dirname(__file__)); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.path.abspath(__file__)).rsplit("/",1)[0]; RES=os.path.join(DATA,"results")
def find(tag,kind):
    hits=glob.glob(os.path.join(DATA,"**",f"*{tag}*{kind}*"),recursive=True)
    if not hits: raise FileNotFoundError(f"missing {tag} {kind}")
    return hits[0]
def load(tag,geno,tp):
    m=sio.mmread(gzip.open(find(tag,"matrix"))).tocsr()
    feats=[l.split("\t")[1] for l in gzip.open(find(tag,"features"),"rt").read().splitlines()]
    bcs=gzip.open(find(tag,"barcodes"),"rt").read().splitlines()
    a=sc.AnnData(m.T.tocsr().astype("float32")); a.var_names=feats; a.var_names_make_unique()
    a.obs_names=[f"{geno}_{tp}_{b}" for b in bcs]; a.obs["genotype"]=geno; a.obs["timepoint"]=tp
    a.var["mt"]=a.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True,percent_top=None)
    a=a[a.obs["pct_counts_mt"]<10].copy(); sc.pp.filter_cells(a,min_genes=300)
    return a
samples=[("wk04_wt","WT","wk04"),("wk04_mut","KO","wk04"),("wk36_wt","WT","wk36"),("wk36_mut","KO","wk36")]
ads=[load(*s) for s in samples]
ad=sc.concat(ads,join="inner")
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
sc.tl.score_genes(ad,[g for g in PROLIFERATION_CORE if g in ad.var_names],score_name="prolif")
# clustering on cell-cycle-regressed HVG (clustering only)
sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor="seurat")
mark={"HSC":["Procr","Hlf","Mecom","Fgd5","Mllt3","Hoxa9","Meis1","Mpl","Ctnnal1"],
 "MkP":["Pf4","Vwf","Itga2b","Gp1bb","Plek"],"Ery":["Klf1","Gata1","Car1","Car2","Gypa"],
 "GMP":["Elane","Mpo","Prtn3","Ctsg","Ms4a3"],"Lymph":["Dntt","Il7r","Vpreb1","Cd79a"]}
for k,gs in mark.items(): sc.tl.score_genes(ad,[g for g in gs if g in ad.var_names],score_name="sc_"+k)
hv=ad.var_names[ad.var.highly_variable]; X=ad[:,hv].X.toarray()
C=np.column_stack([np.ones(ad.n_obs),ad.obs["S_score"],ad.obs["G2M_score"]])
Xr=X-C@np.linalg.lstsq(C,X,rcond=None)[0]
Xs=np.clip((Xr-Xr.mean(0))/(Xr.std(0)+1e-9),-10,10)
pcs=PCA(30,random_state=0).fit_transform(Xs)
ad.obs["cluster"]=KMeans(12,random_state=0,n_init=10).fit_predict(pcs).astype(str)
sig=ad.obs.groupby("cluster")[["sc_"+k for k in mark]].mean()
sigz=(sig-sig.mean())/sig.std(); lab=sigz.idxmax(1).str.replace("sc_","").to_dict()
ad.obs["celltype"]=ad.obs["cluster"].map(lab); ad.obs["cycling"]=ad.obs["phase"]!="G1"
prim=[c for c in sigz.index if sigz.loc[c].idxmax()=="sc_HSC"]
ad.obs["compartment"]=np.where(ad.obs["cluster"].isin(prim),"LT-HSC+MPP1","Committed")
ad.obs.to_csv(os.path.join(RES,"asxl1_10x_obs.csv"))
P=ad.obs[ad.obs.compartment=="LT-HSC+MPP1"]
print("LT-HSC+MPP1 cells per genotype x timepoint:"); print(P.groupby(["timepoint","genotype"]).size().to_string())
# within-compartment, per timepoint (batch-free: WT & KO co-captured) + pooled
def comp(d,lab):
    w=d[d.genotype=="WT"]; k=d[d.genotype=="KO"]
    if len(w)<5 or len(k)<5: return
    chi2,p,_,_=stats.chi2_contingency([[w.cycling.sum(),(~w.cycling).sum()],[k.cycling.sum(),(~k.cycling).sum()]])
    u,pu=stats.mannwhitneyu(k.prolif,w.prolif)
    print(f"  [{lab}] cyc WT={w.cycling.mean()*100:.1f}% KO={k.cycling.mean()*100:.1f}% chi2p={p:.2g} | prolif WT={w.prolif.mean():.3f} KO={k.prolif.mean():.3f} MWUp={pu:.2g}")
print("\nWithin LT-HSC+MPP1 (cell-level):")
for tp in ["wk04","wk36"]: comp(P[P.timepoint==tp],tp)
comp(P,"pooled")
# replicate-aware: pseudobulk per (genotype x timepoint), paired across timepoints
pbk=P.groupby(["genotype","timepoint"]).agg(cyc=("cycling","mean"),prolif=("prolif","mean")).reset_index()
print("\nPseudobulk per genotype x timepoint (replicate-aware units):"); print(pbk.round(3).to_string(index=False))
piv=pbk.pivot(index="timepoint",columns="genotype",values="cyc")
if piv.shape==(2,2):
    t,p=stats.ttest_rel(piv["KO"],piv["WT"]); print(f"\nPaired (timepoint) test KO vs WT cycling: dWT->KO mean={ (piv['KO']-piv['WT']).mean():+.3f}  paired-p={p:.3g} (n=2 pairs)")
print("\nDONE")
