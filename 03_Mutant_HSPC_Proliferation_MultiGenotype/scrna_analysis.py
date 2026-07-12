import os, gzip, tarfile, sys, warnings, io
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats, io as sio
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
from genesets_mouse import S_PHASE, G2M_PHASE, PROLIFERATION_CORE
sc.settings.verbosity = 0
DATA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(DATA, "results"); os.makedirs(RES, exist_ok=True)

def score_and_phase(ad):
    sc.pp.filter_cells(ad, min_genes=300)
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    s = [g for g in S_PHASE if g in ad.var_names]
    g2 = [g for g in G2M_PHASE if g in ad.var_names]
    sc.tl.score_genes_cell_cycle(ad, s_genes=s, g2m_genes=g2)
    pc = [g for g in PROLIFERATION_CORE if g in ad.var_names]
    sc.tl.score_genes(ad, pc, score_name="prolif_core")
    return ad

def summarize(ad, label):
    ph = ad.obs["phase"].value_counts(normalize=True)*100
    cyc = (ad.obs["phase"]!="G1").mean()*100
    return dict(group=label, n_cells=ad.n_obs,
                pct_G1=ph.get("G1",0), pct_S=ph.get("S",0), pct_G2M=ph.get("G2M",0),
                pct_cycling_S_G2M=cyc, mean_prolif_core=ad.obs["prolif_core"].mean())

rows=[]; obs_store={}

# ---------------- GSE124822 Smart-seq2 (WT / Tet2 / Dnmt3a) ----------------
plates = {"WT":["GSE124822_LTHSC-WT_476_RNA_counts.txt.gz","GSE124822_LTHSC-WT_846_RNA_counts.txt.gz"],
          "Tet2":["GSE124822_LTHSC-TET2_30_RNA_counts.txt.gz","GSE124822_LTHSC_TET2_272_RNA_counts.txt.gz"],
          "Dnmt3a":["GSE124822_LTHSC-DNMT3A_675_RNA_counts.txt.gz"]}
for geno, files in plates.items():
    mats=[]
    for f in files:
        d = pd.read_csv(os.path.join(DATA,f), sep="\t", index_col=0)
        d = d.groupby(d.index).sum()
        mats.append(d)
    genes = sorted(set().union(*[set(m.index) for m in mats]))
    mats = [m.reindex(genes).fillna(0) for m in mats]
    expr = pd.concat(mats, axis=1)
    ad = sc.AnnData(expr.T.astype("float32")); ad.var_names_make_unique()
    ad = score_and_phase(ad)
    obs_store[("smartseq2",geno)] = ad.obs.copy()
    rows.append({"dataset":"GSE124822_SmartSeq2", **summarize(ad, geno)})

# ---------------- GSE263101 10x (WT vs Jak2V617F) ----------------
def load_10x(prefix):
    with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
        mtx = sio.mmread(gzip.open(tf.extractfile(prefix+"matrix.mtx.gz"))).tocsr()  # genes x cells
        feats=[l.split("\t")[1] for l in gzip.open(tf.extractfile(prefix+"features.tsv.gz"),"rt").read().splitlines()]
        bcs=[l for l in gzip.open(tf.extractfile(prefix+"barcodes.tsv.gz"),"rt").read().splitlines()]
    ad = sc.AnnData(mtx.T.tocsr().astype("float32")); ad.var_names=feats; ad.obs_names=bcs
    ad.var_names_make_unique()
    return ad
for geno,pref in {"WT":"GSM8185789_Vav1_iCre_","Jak2V617F":"GSM8185791_Jak2V617F_"}.items():
    ad = load_10x(pref)
    ad.var["mt"]=ad.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    ad = ad[ad.obs["pct_counts_mt"]<10].copy()
    ad = score_and_phase(ad)
    obs_store[("tenx",geno)] = ad.obs.copy()
    rows.append({"dataset":"GSE263101_10x", **summarize(ad, geno)})

summ = pd.DataFrame(rows)
summ.to_csv(os.path.join(RES,"scrna_cellcycle_summary.csv"), index=False)
pd.set_option("display.width",170,"display.max_columns",20)
print(summ.to_string(index=False))

# ---- stats: cycling fraction mutant vs WT (chi-square) + prolif score (MWU) ----
def comp(ds, wt, mut, store_keys):
    wo=obs_store[store_keys[0]]; mo=obs_store[store_keys[1]]
    cyc_w=(wo["phase"]!="G1"); cyc_m=(mo["phase"]!="G1")
    tab=np.array([[cyc_m.sum(), (~cyc_m).sum()],[cyc_w.sum(), (~cyc_w).sum()]])
    chi2,p,_,_=stats.chi2_contingency(tab)
    u,pu=stats.mannwhitneyu(mo["prolif_core"], wo["prolif_core"], alternative="two-sided")
    return dict(dataset=ds, comparison=f"{mut} vs {wt}",
                cycling_WT=cyc_w.mean()*100, cycling_MUT=cyc_m.mean()*100,
                chi2_p=p, prolif_WT=wo["prolif_core"].mean(), prolif_MUT=mo["prolif_core"].mean(), MWU_p=pu)
stats_rows=[comp("GSE124822","WT","Tet2",[("smartseq2","WT"),("smartseq2","Tet2")]),
            comp("GSE124822","WT","Dnmt3a",[("smartseq2","WT"),("smartseq2","Dnmt3a")]),
            comp("GSE263101","WT","Jak2V617F",[("tenx","WT"),("tenx","Jak2V617F")])]
sdf=pd.DataFrame(stats_rows); sdf.to_csv(os.path.join(RES,"scrna_comparisons.csv"),index=False)
print("\n=== comparisons ==="); print(sdf.to_string(index=False))
