#!/usr/bin/env python3
# CH / myeloid-inflammation gene heatmap across the 4 BMDM groups.
#   Panel A: expression level (z-scored mean log2 CPM) per group.
#   Panel B: signed effect (mean log2-CPM difference) for the 4 pairwise
#            contrasts, with BH-corrected significance stars.
#     burden (HB vs LB) -> UNPAIRED t-test
#     lenti  (pos vs neg, donor-matched) -> PAIRED t-test
# Each test uses the CPM file where both its groups were co-normalized.

import os, re
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize

BMDM = "/sessions/dreamy-exciting-darwin/mnt/BMDM"
OUT  = "/sessions/dreamy-exciting-darwin/mnt/outputs"

# ---- gene panel (grouped) ----
GENES = (["Il1b","Il1a","Il6","Tnf","Il18"]                      # cytokines
       + ["Nlrp3","Casp1","Pycard","Aim2"]                       # inflammasome
       + ["Ccl2","Ccl5","Cxcl1","Cxcl2","Cxcl10"]                # chemokines
       + ["Tlr2","Tlr4","Nfkb1","Tnfaip3","Ptgs2","Nos2"])       # signaling/effectors

COMPS = {
    "HB.neg_vs_LB.neg": dict(pos=None),
    "HB.pos_vs_HB.neg": dict(),
    "HB.pos_vs_LB.pos": dict(),
    "LB.pos_vs_LB.neg": dict(),
}

def deg_path(f):  return os.path.join(BMDM, f, f"DEG_{f}_bulk.txt")
def cpm_path(f):  return os.path.join(BMDM, f, f"edgeR_{f}_TMM_normalized_CPM_wOffset.csv")

# ---- symbol -> ensembl map (union across DEG files) ----
sym2ens = {}
for f in COMPS:
    d = pd.read_csv(deg_path(f), sep="\t")
    d = d.rename(columns={d.columns[0]: "ensembl"})
    for ens, sym in zip(d["ensembl"], d["gene_name"]):
        sym2ens.setdefault(sym, ens)
missing = [g for g in GENES if g not in sym2ens]
if missing: print("WARNING: not in any DEG file:", missing)

def load_cpm(f):
    c = pd.read_csv(cpm_path(f), index_col=0)
    return c  # rows=ensembl, cols=samples

cpm = {f: load_cpm(f) for f in COMPS}

def cols(df, pat):       return [c for c in df.columns if re.search(pat, c)]
def log2cpm(df, samples, ens):
    if ens not in df.index: return None
    return np.log2(df.loc[ens, samples].astype(float).values + 1.0)

# group-level mean log2 CPM (neg groups from file1, pos groups from file3)
f_neg = "HB.neg_vs_LB.neg"; f_pos = "HB.pos_vs_LB.pos"
group_src = {
    "HB.lentineg": (f_neg, r"^HB\d+neg$"),
    "HB.lentipos": (f_pos, r"^HB\d+pos$"),
    "LB.lentineg": (f_neg, r"^LB\d+neg$"),
    "LB.lentipos": (f_pos, r"^LB\d+pos$"),
}
GROUP_ORDER = ["HB.lentineg","HB.lentipos","LB.lentineg","LB.lentipos"]

mean_log2, mean_cpm = {}, {}
for g in GROUP_ORDER:
    f, pat = group_src[g]
    df = cpm[f]; sc = cols(df, pat)
    mean_log2[g] = {}; mean_cpm[g] = {}
    for sym in GENES:
        ens = sym2ens.get(sym)
        v = log2cpm(df, sc, ens) if ens else None
        mean_log2[g][sym] = np.nan if v is None else v.mean()
        mean_cpm[g][sym]  = np.nan if (ens is None or ens not in df.index) else df.loc[ens, sc].astype(float).mean()

L = pd.DataFrame(mean_log2)[GROUP_ORDER].reindex(GENES)     # mean log2 CPM
C = pd.DataFrame(mean_cpm)[GROUP_ORDER].reindex(GENES)      # mean linear CPM
# z-score per gene across the 4 groups
Z = L.sub(L.mean(axis=1), axis=0).div(L.std(axis=1, ddof=0).replace(0, np.nan), axis=0)

# ---- pairwise contrasts ----
# (label, file, groupA_pat, groupB_pat, paired, effect = A - B; sign convention)
CONTRASTS = [
    ("HB vs LB | neg\n(unpaired)",  f_neg,              r"^HB\d+neg$", r"^LB\d+neg$", False),
    ("HB vs LB | pos\n(unpaired)",  f_pos,              r"^HB\d+pos$", r"^LB\d+pos$", False),
    ("pos vs neg | HB\n(paired)",   "HB.pos_vs_HB.neg", r"^HB\d+pos$", r"^HB\d+neg$", True),
    ("pos vs neg | LB\n(paired)",   "LB.pos_vs_LB.neg", r"^LB\d+pos$", r"^LB\d+neg$", True),
]
def donor(s): return re.sub(r"(pos|neg)$", "", s)

eff = pd.DataFrame(index=GENES, columns=[c[0] for c in CONTRASTS], dtype=float)
pval = eff.copy()
for label, f, pa, pb, paired in CONTRASTS:
    df = cpm[f]; a = cols(df, pa); b = cols(df, pb)
    if paired:  # align donors
        da = {donor(x): x for x in a}; db = {donor(x): x for x in b}
        common = [d for d in da if d in db]
        a = [da[d] for d in common]; b = [db[d] for d in common]
    for sym in GENES:
        ens = sym2ens.get(sym)
        va = log2cpm(df, a, ens); vb = log2cpm(df, b, ens)
        if va is None or vb is None: continue
        eff.loc[sym, label] = va.mean() - vb.mean()
        if paired:
            t, p = stats.ttest_rel(va, vb)
        else:
            t, p = stats.ttest_ind(va, vb, equal_var=False)  # Welch
        pval.loc[sym, label] = p

# BH-FDR within each contrast (20 genes)
def bh(p):
    p = np.asarray(p, float); ok = ~np.isnan(p); q = np.full_like(p, np.nan)
    pp = p[ok]; n = len(pp); order = np.argsort(pp); ranked = pp[order]
    qv = ranked * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(qv, 0, 1); q[ok] = out
    return q
qval = pd.DataFrame({c: bh(pval[c].values) for c in pval.columns}, index=GENES)

def star(q):
    if np.isnan(q): return ""
    return "***" if q<.001 else "**" if q<.01 else "*" if q<.05 else "·" if q<.10 else ""

# ---- save stats table ----
stats_out = pd.DataFrame(index=GENES)
stats_out["ensembl"] = [sym2ens.get(g) for g in GENES]
for g in GROUP_ORDER:
    stats_out[f"meanCPM_{g}"] = C[g]
for label, *_ in [(c[0],) for c in CONTRASTS]:
    short = label.split("\n")[0].replace(" ","").replace("|","_")
    stats_out[f"log2diff_{short}"] = eff[label]
    stats_out[f"p_{short}"]  = pval[label]
    stats_out[f"BHq_{short}"] = qval[label]
stats_out.to_csv(os.path.join(OUT, "ch_myeloid_heatmap_stats.csv"))

# ============================ FIGURE ============================
ng = len(GENES)
fig = plt.figure(figsize=(13.5, 0.46*ng + 2.6))
gs = fig.add_gridspec(1, 2, width_ratios=[4, 4], wspace=0.55)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])

# group block boundaries (cytokine/inflammasome/chemokine/signaling)
block_after = [5, 9, 14]  # draw lines after these row indices

# ---- Panel A : expression level (z-scored) ----
zmax = np.nanmax(np.abs(Z.values))
imA = axA.imshow(Z.values, aspect="auto", cmap="RdBu_r",
                 norm=TwoSlopeNorm(vmin=-zmax, vcenter=0, vmax=zmax))
axA.set_xticks(range(4)); axA.set_xticklabels(GROUP_ORDER, rotation=30, ha="right", fontsize=9)
axA.set_yticks(range(ng)); axA.set_yticklabels(GENES, fontstyle="italic", fontsize=9)
axA.set_title("Expression level\n(row z-score of mean log2 CPM)", fontsize=10)
for i in range(ng):
    for j in range(4):
        v = C.values[i, j]
        if not np.isnan(v):
            axA.text(j, i, f"{v:.0f}" if v>=10 else f"{v:.1f}",
                     ha="center", va="center", fontsize=6.5,
                     color="black" if abs(Z.values[i,j])<1.3 else "white")
for b in block_after:
    axA.axhline(b-0.5, color="0.25", lw=0.8)
axA.set_xticks(np.arange(-.5,4,1), minor=True); axA.set_yticks(np.arange(-.5,ng,1), minor=True)
axA.grid(which="minor", color="white", lw=0.6); axA.tick_params(which="minor", length=0)
cbA = fig.colorbar(imA, ax=axA, fraction=0.046, pad=0.02); cbA.set_label("z-score", fontsize=8)

# ---- Panel B : effect (log2 diff) + stars ----
emax = np.nanmax(np.abs(eff.values))
imB = axB.imshow(eff.values, aspect="auto", cmap="PRGn",
                 norm=TwoSlopeNorm(vmin=-emax, vcenter=0, vmax=emax))
clabels = [c[0] for c in CONTRASTS]
axB.set_xticks(range(4)); axB.set_xticklabels(clabels, fontsize=8)
axB.set_yticks(range(ng)); axB.set_yticklabels([])
axB.set_title("Effect size & significance\nlog2 mean-CPM difference (A−B), BH-FDR stars", fontsize=10)
for i in range(ng):
    for j in range(4):
        v = eff.values[i, j]; q = qval.values[i, j]
        if np.isnan(v): continue
        s = star(q)
        axB.text(j, i, s, ha="center", va="center", fontsize=11, fontweight="bold",
                 color="white" if abs(v) > 0.5*emax else "black")
axB.axvline(1.5, color="0.25", lw=1.0)  # split burden | lenti
for b in block_after:
    axB.axhline(b-0.5, color="0.25", lw=0.8)
axB.set_xticks(np.arange(-.5,4,1), minor=True); axB.set_yticks(np.arange(-.5,ng,1), minor=True)
axB.grid(which="minor", color="white", lw=0.6); axB.tick_params(which="minor", length=0)
cbB = fig.colorbar(imB, ax=axB, fraction=0.046, pad=0.02)
cbB.set_label("log2 difference\n(+ = up in HB / up in pos)", fontsize=8)

fig.suptitle("Clonal-hematopoiesis myeloid-inflammation genes across BMDM groups",
             fontsize=13, y=0.995)
fig.text(0.5, 0.005,
         "Stars: *** BH-q<0.001  ** <0.01  * <0.05  · <0.10   |   "
         "burden = unpaired t-test, lenti = donor-paired t-test (n: HB=5, LB=4)   |   "
         "neg groups co-normalized in HB.neg_vs_LB.neg; pos groups in HB.pos_vs_LB.pos",
         ha="center", fontsize=7.2, color="0.30")
fig.savefig(os.path.join(OUT, "ch_myeloid_heatmap.png"), dpi=170, bbox_inches="tight")
print("Wrote ch_myeloid_heatmap.png and ch_myeloid_heatmap_stats.csv")

# quick console summary of significant hits
for label, *_ in [(c[0],) for c in CONTRASTS]:
    sig = [(g, eff.loc[g,label], qval.loc[g,label]) for g in GENES
           if not np.isnan(qval.loc[g,label]) and qval.loc[g,label] < 0.05]
    print(f"\n{label.split(chr(10))[0]}  -> {len(sig)} genes BH-q<0.05")
    for g,e,q in sorted(sig, key=lambda x:x[2]):
        print(f"   {g:8s} log2diff={e:+.2f}  q={q:.3g}")
