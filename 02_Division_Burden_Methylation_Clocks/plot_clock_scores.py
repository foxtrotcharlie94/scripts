import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('/tmp/gse44117/clock_cpg_analysis.csv')

group_order = ["young_baseline", "old_baseline", "young_10_reconst", "old_10_reconst"]
group_labels = {
    "young_baseline": "Young HSC\n(baseline, n=2)",
    "old_baseline": "Old HSC\n(baseline, n=2)",
    "young_10_reconst": "Young donor,\n10-cell transplant\n(n=3)",
    "old_10_reconst": "Old-recipient,\n10-cell transplant\n(n=2)\n[not in Fig 2E]",
}
colors = {"young_baseline": "#4C72B0", "old_baseline": "#DD8452",
          "young_10_reconst": "#55A868", "old_10_reconst": "#C44E52"}

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for ax, clockname in zip(axes, ["Blood", "WLMT"]):
    sub = df[df["clock"] == clockname]
    grp = sub.groupby("group")["weighted_sum_score_normalized"].agg(["mean", "std", "count"]).reindex(group_order)
    x = np.arange(len(grp))
    ax.bar(x, grp["mean"], yerr=grp["std"], capsize=5,
           color=[colors[g] for g in grp.index], alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([group_labels[g] for g in grp.index], fontsize=8)
    ax.axhline(0, color='gray', lw=0.8)
    ax.set_ylabel("Weighted clock score\n(weight-normalized sum over covered CpGs)")
    ax.set_title(f"{clockname} clock\n(CpG sites lifted mm10→mm9, matched to RRBS)")

plt.suptitle("GSE44117 (Beerman et al. 2013) HSC RRBS\nActual clock-CpG weighted scores by group", y=1.02)
plt.tight_layout()
plt.savefig('/tmp/gse44117/clock_scores_by_group.png', dpi=150, bbox_inches='tight')
print("saved plot")

print(df.groupby(["clock","group"])["weighted_sum_score_normalized"].agg(["mean","std","count"]).reindex(
    pd.MultiIndex.from_product([["Blood","WLMT"], group_order])
))

# Also a positive vs negative weight cpg view
fig2, axes2 = plt.subplots(1, 2, figsize=(11, 5))
for ax, clockname in zip(axes2, ["Blood", "WLMT"]):
    sub = df[df["clock"] == clockname]
    grp = sub.groupby("group")[["mean_pct_meth_pos_weight_cpgs","mean_pct_meth_neg_weight_cpgs"]].mean().reindex(group_order)
    x = np.arange(len(grp))
    width = 0.35
    ax.bar(x - width/2, grp["mean_pct_meth_pos_weight_cpgs"], width, label="Positive-weight CpGs\n(meth ↑ with age)", color="#2ca02c", alpha=0.85)
    ax.bar(x + width/2, grp["mean_pct_meth_neg_weight_cpgs"], width, label="Negative-weight CpGs\n(meth ↓ with age)", color="#d62728", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([group_labels[g] for g in grp.index], fontsize=8)
    ax.set_ylabel("Mean % methylation")
    ax.set_title(f"{clockname} clock CpGs")
    ax.legend(fontsize=7)
plt.suptitle("Mean methylation at clock CpGs, split by weight sign", y=1.02)
plt.tight_layout()
plt.savefig('/tmp/gse44117/clock_cpg_posneg_by_group.png', dpi=150, bbox_inches='tight')
print("saved plot 2")
