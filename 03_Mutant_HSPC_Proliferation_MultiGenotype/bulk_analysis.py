import os, io, gzip, tarfile, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from genesets_mouse import GENESETS

DATA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(DATA, "results")
os.makedirs(RES, exist_ok=True)

def cpm_log(df):
    cpm = df / df.sum(axis=0) * 1e6
    return np.log2(cpm + 1)

def zscore_rows(df):
    m = df.mean(axis=1).values[:,None]; s = df.std(axis=1, ddof=0).values[:,None]
    s[s==0] = 1.0
    return pd.DataFrame((df.values - m)/s, index=df.index, columns=df.columns)

def sig_scores(z, genes):
    idx = {g.lower(): g for g in z.index}
    use = [idx[g.lower()] for g in genes if g.lower() in idx]
    return z.loc[use].mean(axis=0), len(use)

def analyze(name, logexpr, wt, mut):
    # collapse duplicate gene rows by mean
    logexpr = logexpr.groupby(logexpr.index).mean()
    z = zscore_rows(logexpr[wt+mut])
    rows = []
    per_sample = {}
    for setname, genes in GENESETS.items():
        sc, n = sig_scores(z, genes)
        per_sample[setname] = sc
        a = sc[mut].values; b = sc[wt].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append(dict(dataset=name, geneset=setname, n_genes=n,
                         WT_mean=b.mean(), MUT_mean=a.mean(),
                         diff_MUT_minus_WT=a.mean()-b.mean(), t=t, p=p))
    ss = pd.DataFrame(per_sample); ss.index.name="sample"
    ss["group"] = ["WT" if s in wt else "MUT" for s in ss.index]
    ss.to_csv(os.path.join(RES, f"{name}_signature_scores_persample.csv"))
    return pd.DataFrame(rows), ss

summary = []

# ---------- TET2  GSE96758 (already normalized, log-scale) ----------
t2 = pd.read_csv(os.path.join(DATA,"TET2GSE96758_normalizedData.txt.gz"),
                 sep="\t", index_col=0)
wt=["WT1","WT2","WT3"]; mut=["T1","T2","T3"]
r,_ = analyze("Tet2_GSE96758", t2, wt, mut); summary.append(r)

# ---------- ASXL1 GSE158184 (raw counts) ----------
ax = pd.read_csv(os.path.join(DATA,"GSE158184_Bulk_ASXL_counts.filt.txt.gz"),
                 sep="\t", index_col=0)
wt=["W2","W3","W4","W5"]; mut=["M2","M3","M4","M5"]
r,_ = analyze("Asxl1_GSE158184", cpm_log(ax[wt+mut]), wt, mut); summary.append(r)

# ---------- JAK2 GSE123401 (featureCounts, Ensembl ids) ----------
# Ensembl->symbol map from the 10x features file in the other Jak2 dataset
feat = {}
with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
    f = tf.extractfile("GSM8185789_Vav1_iCre_features.tsv.gz")
    for line in gzip.open(f, "rt"):
        e,s = line.rstrip("\n").split("\t")[:2]; feat[e]=s
jmap = {"GSM3502680":"WT1","GSM3502681":"WT2","GSM3502682":"WT3",
        "GSM3502683":"VF1","GSM3502684":"VF2","GSM3502685":"VF3"}
cols = {}
with tarfile.open(os.path.join(DATA,"GSE123401_RAW.tar")) as tf:
    for m in tf.getmembers():
        gsm = m.name.split("_")[0]
        if gsm in jmap:
            df = pd.read_csv(gzip.open(tf.extractfile(m),"rt"), sep="\t", comment="#")
            cnt = df.iloc[:,[0,-1]]; cnt.columns=["gid","count"]
            cnt["sym"] = cnt["gid"].str.split(".").str[0].map(feat)
            cnt = cnt.dropna(subset=["sym"]).groupby("sym")["count"].sum()
            cols[jmap[gsm]] = cnt
jak = pd.DataFrame(cols)
wt=["WT1","WT2","WT3"]; mut=["VF1","VF2","VF3"]
r,_ = analyze("Jak2V617F_GSE123401", cpm_log(jak[wt+mut]), wt, mut); summary.append(r)

out = pd.concat(summary, ignore_index=True)
out.to_csv(os.path.join(RES,"bulk_proliferation_summary.csv"), index=False)
pd.set_option("display.width",160, "display.max_columns",20)
print(out.to_string(index=False))
