import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import re

df = pd.read_csv("/tmp/gse44117/summary.csv")

def group_of(sample):
    if sample.startswith("HSC_fetal"):
        return "fetal"
    m = re.match(r"HSC_(young|old)_(\d+)_reconst", sample)
    if m:
        return f"{m.group(1)}_{m.group(2)}_reconst"
    m = re.match(r"HSC_(young)_(\d)x5FU", sample)
    if m:
        return f"young_{m.group(2)}x5FU"
    m = re.match(r"HSC_(young|old)_\d+$", sample)
    if m:
        return f"{m.group(1)}_baseline"
    return "other"

df["group"] = df["sample"].apply(group_of)
df.to_csv("/tmp/gse44117/summary_with_groups.csv", index=False)

# Group-level summary (mean +/- SD of weighted % methylation, n_cpgs, coverage)
grp = df.groupby("group").agg(
    n_samples=("sample", "count"),
    mean_pct_meth=("mean_pct_meth_weighted", "mean"),
    sd_pct_meth=("mean_pct_meth_weighted", "std"),
    mean_n_cpgs=("n_cpgs", "mean"),
    mean_coverage=("mean_coverage", "mean"),
).reset_index()
grp.to_csv("/tmp/gse44117/group_summary.csv", index=False)
print(grp.to_string(index=False))

# Figure 2E-relevant groups only: baseline young/old + 10-cell reconstitution young/old
fig2e_groups = ["young_baseline", "old_baseline", "young_10_reconst", "old_10_reconst"]
sub = grp[grp["group"].isin(fig2e_groups)].copy()
order = {g: i for i, g in enumerate(fig2e_groups)}
sub["order"] = sub["group"].map(order)
sub = sub.sort_values("order")
print("\nFig 2E-relevant subset:")
print(sub.to_string(index=False))

labels = {
    "young_baseline": "Young HSC\n(baseline, n=2)",
    "old_baseline": "Old HSC\n(baseline, n=2)",
    "young_10_reconst": "Young donor\n10-cell transplant\n(20wk engraft, n=3)",
    "old_10_reconst": "Old recipient\n10-cell transplant\n(20wk engraft, n=2)",
}

fig, ax = plt.subplots(figsize=(6, 5))
x = np.arange(len(sub))
ax.bar(x, sub["mean_pct_meth"], yerr=sub["sd_pct_meth"], capsize=5,
       color=["#4C72B0", "#DD8452", "#4C72B0", "#DD8452"], alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels([labels[g] for g in sub["group"]], fontsize=9)
ax.set_ylabel("Global mean CpG methylation (%)\n(weighted, RRBS)")
ax.set_title("GSE44117 (Beerman et al. 2013) HSC RRBS\nGlobal methylation by group\n(proxy only — NOT the BS-WLMT/BS-Mouse Blood clock scores in Fig 2E)")
plt.tight_layout()
plt.savefig("/tmp/gse44117/fig2e_groups_global_methylation.png", dpi=150)
print("\nSaved plot.")
