#!/usr/bin/env python3
# BMDM asymmetry diagnostics (works from CPM + DEG only; no raw counts).
#   1) CPM-based library-composition check: how much of each sample's
#      (TMM-normalized) library is eaten by its top genes.  CPM/1e6 = share.
#   2) logFC distribution + tail asymmetry per comparison.
# Outputs two PNGs + a composition CSV.

import os, re, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D   = "/sessions/dreamy-exciting-darwin/mnt/BMDM"
OUT = "/sessions/dreamy-exciting-darwin/mnt/outputs"

COMPS = ["HB.neg_vs_LB.neg","HB.pos_vs_LB.pos","HB.pos_vs_HB.neg","LB.pos_vs_LB.neg"]
POS   = {"HB.neg_vs_LB.neg":"+ = up in HB","HB.pos_vs_LB.pos":"+ = up in HB",
         "HB.pos_vs_HB.neg":"+ = up in pos","LB.pos_vs_LB.neg":"+ = up in pos"}

def cpm_path(f): return os.path.join(D, f, f"edgeR_{f}_TMM_normalized_CPM_wOffset.csv")
def deg_path(f): return os.path.join(D, f, f"DEG_{f}_bulk.txt")

# symbol map for labeling top genes
sym = {}
for f in COMPS:
    d = pd.read_csv(deg_path(f), sep="\t").rename(columns={None:"ens"})
    d = d.rename(columns={d.columns[0]:"ens"})
    for e,s in zip(d["ens"], d["gene_name"]): sym.setdefault(e,s)

# ============================================================
# 1) CPM-based composition check
# ============================================================
# Build one column per *sample* (each sample taken from one file).
# neg samples from HB.neg_vs_LB.neg ; pos samples from HB.pos_vs_LB.pos
cN = pd.read_csv(cpm_path("HB.neg_vs_LB.neg"), index_col=0)
cP = pd.read_csv(cpm_path("HB.pos_vs_LB.pos"), index_col=0)
samples = {}
for s in cN.columns: samples[s] = cN[s]
for s in cP.columns: samples[s] = cP[s]
M = pd.DataFrame(samples)  # genes x samples (CPM)

def group_of(s):
    cond = "HB" if s.startswith("HB") else "LB"
    geno = "pos" if s.endswith("pos") else "neg"
    return f"{cond}.lenti{geno}"
GROUPS = ["HB.lentineg","HB.lentipos","LB.lentineg","LB.lentipos"]
samp_by_group = {g:[s for s in M.columns if group_of(s)==g] for g in GROUPS}

rows = []
for s in M.columns:
    v = M[s].sort_values(ascending=False)
    share = v / 1e6                      # fraction of library
    rows.append(dict(sample=s, group=group_of(s),
        top1=share.iloc[0], top1_gene=sym.get(v.index[0], v.index[0]),
        top5=share.iloc[:5].sum(), top10=share.iloc[:10].sum(),
        top20=share.iloc[:20].sum(), top50=share.iloc[:50].sum()))
comp = pd.DataFrame(rows)
comp.to_csv(os.path.join(OUT, "bmdm_composition_check.csv"), index=False)

print("=== Library composition (fraction of TMM-normalized CPM in top-N genes) ===")
print("Group means:")
print((comp.groupby("group")[["top1","top5","top10","top20","top50"]].mean()
       .reindex(GROUPS).round(3)).to_string())
print("\nMost dominant single gene per sample:")
print(comp[["sample","top1","top1_gene"]].sort_values("top1",ascending=False).head(8).to_string(index=False))

# group-mean top-gene table (which genes dominate, by group)
print("\nTop 12 genes by mean CPM share, per group:")
group_top = {}
for g in GROUPS:
    gm = M[samp_by_group[g]].mean(axis=1).sort_values(ascending=False)
    group_top[g] = [(sym.get(i,i), gm[i]/1e6) for i in gm.index[:12]]
for g in GROUPS:
    print(f"  {g}: " + ", ".join(f"{n}({s*100:.1f}%)" for n,s in group_top[g][:8]))

# ---- composition figure: stacked top-N share per sample ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axes[0]
order = [s for g in GROUPS for s in samp_by_group[g]]
gcol = {"HB.lentineg":"#4C72B0","HB.lentipos":"#8FB0DD",
        "LB.lentineg":"#C44E52","LB.lentipos":"#E69598"}
xs = np.arange(len(order))
for k,(lvl,lab) in enumerate([("top1","top 1"),("top5","top 5"),("top10","top 10"),
                              ("top20","top 20"),("top50","top 50")]):
    vals = comp.set_index("sample").loc[order, lvl].values*100
    ax.plot(xs, vals, marker="o", label=lab)
ax.set_xticks(xs); ax.set_xticklabels(order, rotation=90, fontsize=7)
ax.set_ylabel("% of library (TMM-normalized CPM)")
ax.set_title("Library dominance by top-N genes\n(per sample)")
ax.legend(fontsize=8, title="cumulative")
ax.grid(axis="y", alpha=0.3)
for g in GROUPS:  # shade group bands
    idx=[order.index(s) for s in samp_by_group[g]]
    ax.axvspan(min(idx)-0.5, max(idx)+0.5, color=gcol[g], alpha=0.06)

# group-mean cumulative-share curve (how fast library accumulates)
ax2 = axes[1]
for g in GROUPS:
    gm = M[samp_by_group[g]].mean(axis=1).sort_values(ascending=False).values/1e6
    cum = np.cumsum(gm)*100
    ax2.plot(np.arange(1,len(cum)+1), cum, label=g, color=gcol[g], lw=1.8)
ax2.set_xscale("log"); ax2.set_xlim(1, 2000)
ax2.set_xlabel("number of top genes (log)"); ax2.set_ylabel("cumulative % of library")
ax2.set_title("Cumulative library share (group mean)")
ax2.axhline(50, ls="--", color="grey", lw=0.8); ax2.text(1.1,51,"50%",fontsize=7,color="grey")
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
fig.suptitle("BMDM library-composition check (from CPM; TMM already applied)", fontsize=12)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(os.path.join(OUT,"bmdm_composition_check.png"), dpi=160, bbox_inches="tight")
print("\nWrote bmdm_composition_check.png / .csv")

# ============================================================
# 2) logFC distribution + tail asymmetry
# ============================================================
fig, axes = plt.subplots(2, 4, figsize=(18, 8))
print("\n=== logFC distribution / tail asymmetry ===")
for j,c in enumerate(COMPS):
    d = pd.read_csv(deg_path(c), sep="\t")
    lfc = d["logFC"].values
    med = np.median(lfc); mean = lfc.mean()
    # top tail by |logFC|
    t = d.reindex(d["logFC"].abs().sort_values(ascending=False).index)
    def tailsplit(n):
        h=t.head(n); return (h["logFC"]>0).sum(), (h["logFC"]<0).sum()
    # significant tail (FDR<0.25) split
    sig=d[d["FDR"]<0.25]
    print(f"{c:18s} median={med:+.3f} mean={mean:+.3f} | "
          f"top50 up/dn={tailsplit(50)} top200={tailsplit(200)} | "
          f"FDR<0.25 up/dn=({(sig['logFC']>0).sum()}/{(sig['logFC']<0).sum()}) [{POS[c]}]")

    # --- top: histogram of logFC ---
    axh = axes[0, j]
    axh.hist(lfc, bins=80, color="#888", alpha=0.85)
    axh.axvline(0, color="k", lw=0.8); axh.axvline(med, color="firebrick", lw=1.2, ls="--")
    axh.set_title(c, fontsize=10)
    axh.set_xlabel("logFC"); axh.set_xlim(-3,3)
    axh.text(0.02,0.95, f"median={med:+.3f}\nmean={mean:+.3f}", transform=axh.transAxes,
             va="top", fontsize=8, color="firebrick")
    if j==0: axh.set_ylabel("# genes")

    # --- bottom: ranked logFC (waterfall) showing tail one-sidedness ---
    axr = axes[1, j]
    s = np.sort(lfc)[::-1]
    axr.plot(np.arange(len(s)), s, color="#333", lw=0.7)
    axr.axhline(0, color="grey", lw=0.7)
    axr.fill_between(np.arange(len(s)), s, 0, where=s>0, color="#C44E52", alpha=0.5)
    axr.fill_between(np.arange(len(s)), s, 0, where=s<0, color="#4C72B0", alpha=0.5)
    axr.set_xlabel("gene rank");
    if j==0: axr.set_ylabel("logFC (sorted)")
    axr.set_ylim(-max(abs(s.min()),s.max())*1.05, max(abs(s.min()),s.max())*1.05)
    axr.set_title(POS[c], fontsize=8)
fig.suptitle("logFC distribution (top) & ranked-logFC waterfall (bottom) — tail asymmetry",
             fontsize=13)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(os.path.join(OUT,"bmdm_logFC_distribution.png"), dpi=150, bbox_inches="tight")
print("\nWrote bmdm_logFC_distribution.png")
