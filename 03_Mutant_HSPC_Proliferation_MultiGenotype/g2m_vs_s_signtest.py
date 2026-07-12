"""
Asxl1-KO vs WT — per-cluster sign test using the cell-cycle score split into its
components: S.Score alone, G2M.Score alone, and (S+G2M) combined.
No raw extraction needed — S.Score and G2M.Score are already columns in the metadata.
Unit = mouse (7 WT vs 7 KO). Run after the workspace restarts.
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = "/sessions/<session>/mnt/mutant_proliefration_analysis"  # or the Windows path
m = pd.read_csv(f"{BASE}/seurat_harmony_output/Asxl1/Asxl1_metadata_clusters.csv") \
      .rename(columns={"Unnamed: 0": "barcode"})
m["cc"] = m["S.Score"] + m["G2M.Score"]

def signtest(res, score):
    col = f"RNA_snn_res.{res}"; rows = []
    for cl, g in m.groupby(col):
        pm = g.groupby(["mouse", "genotype"])[score].mean().reset_index()
        wt = pm[pm.genotype == "WT"][score]; ko = pm[pm.genotype == "KO"][score]
        if len(wt) < 2 or len(ko) < 2: continue
        rows.append(ko.mean() - wt.mean())
    d = np.array(rows); up = (d > 0).sum(); tot = len(d)
    sign = stats.binomtest(up, tot, 0.5, alternative="greater").pvalue
    wil = stats.wilcoxon(d).pvalue if tot >= 6 else np.nan
    return up, tot, sign, wil, np.median(d)

print(f"{'metric':16s}{'res':>5}{'up/tot':>9}{'sign p':>9}{'wilcox':>9}{'medianD':>10}")
for res in ["0.4", "1.2"]:
    for s, lab in [("S.Score", "S only"), ("G2M.Score", "G2M only"), ("cc", "S+G2M")]:
        up, tot, sg, wl, md = signtest(res, s)
        print(f"{lab:16s}{res:>5}{f'{up}/{tot}':>9}{sg:>9.3f}{wl:>9.3f}{md:>+10.4f}")
