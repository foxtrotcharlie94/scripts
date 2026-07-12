import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

df1 = pd.read_csv('/tmp/gse44117/cpg_deltas.csv')  # Blood, WLMT
df2 = pd.read_csv('/tmp/gse44117/cpg_deltas_clocks2.csv')
df2 = df2[df2['clock'] == 'YOMT']  # drop Liver
df = pd.concat([df1, df2], ignore_index=True)
df.to_csv('/tmp/gse44117/cpg_deltas_all3.csv', index=False)

fig, ax = plt.subplots(figsize=(8.5, 7.5))
colors = {"Blood": "#4C72B0", "WLMT": "#DD8452", "YOMT": "#55A868"}
text_y = {"Blood": 0.97, "WLMT": 0.90, "YOMT": 0.83}

for clockname in ["Blood", "WLMT", "YOMT"]:
    sub = df[df["clock"] == clockname]
    x, y = sub["delta_age_old_minus_young"], sub["delta_reconst_minus_young"]
    r, p = stats.pearsonr(x, y)
    slope, intercept, r_lin, p_lin, se = stats.linregress(x, y)

    ax.scatter(x, y, s=16, alpha=0.5, color=colors[clockname], label=f"{clockname} clock CpGs (n={len(sub)})")
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, slope * xs + intercept, color=colors[clockname], lw=2)

    ax.text(0.03, text_y[clockname],
            f"{clockname}: slope={slope:.2f}, r={r:.2f}, p={p:.1e}, n={len(sub)}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor=colors[clockname], alpha=0.9))
    print(f"{clockname}: slope={slope:.3f}, r={r:.3f}, p={p:.2e}, n={len(sub)}")

ax.axhline(0, color='lightgray', lw=0.7, zorder=0)
ax.axvline(0, color='lightgray', lw=0.7, zorder=0)
ax.set_xlabel("Δ methylation: Old baseline − Young baseline\n(% points, per CpG, replicate means)")
ax.set_ylabel("Δ methylation: Young 10-cell reconst. − Young baseline\n(% points, per CpG, replicate means)")
ax.set_title("Per-CpG change: natural aging vs. forced-proliferation transplant\n(GSE44117; Blood, WLMT, and YOMT clock CpGs)")
ax.legend(fontsize=8, loc="lower right")
plt.tight_layout()
plt.savefig('/tmp/gse44117/delta_age_vs_delta_reconst_overlay3.png', dpi=150, bbox_inches='tight')
print("saved overlay plot with 3 clocks")
