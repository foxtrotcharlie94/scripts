# Asxl1 reconciliation: is the proliferation signal S-phase or G2M-driven?

**Script:** `g2m_vs_s_signtest.py`

## What it does
A later reconciliation step against an R/Seurat-based Asxl1-KO vs WT clustering (not this project's own Python clustering): reads per-cell Seurat metadata that already contains S.Score/G2M.Score, and for each cluster (at two clustering resolutions, 0.4 and 1.2) computes the mouse-level (7 WT vs 7 KO) mean difference in S-phase score alone, G2M score alone, and the combined S+G2M score, then runs a cross-cluster sign test and Wilcoxon signed-rank test on each of the three metrics separately.

## Inputs
Asxl1_metadata_clusters.csv — Seurat per-cell metadata with columns RNA_snn_res.0.4 / RNA_snn_res.1.2 (cluster IDs), S.Score, G2M.Score, mouse, genotype.

## Outputs
Console table of up/total clusters, sign-test p, Wilcoxon p and median delta for each of S-only / G2M-only / S+G2M, at both resolutions.

## Conceptual aim / target
Determine whether an apparent increase in the combined cell-cycle score in Asxl1-KO reflects a genuine increase in proliferation (both S and G2M phases up) or is driven asymmetrically by just one phase, which would complicate a simple "increased proliferation" interpretation.

## Conclusions / findings
Directly resolves the ambiguity flagged in this project's session history: splitting combined cell-cycle score into components showed the "up" signal in Asxl1-KO progenitors was not a genuine, symmetric proliferation increase — the S-phase-only component behaved differently from G2M, indicating the naive combined-score result should not be read as straightforward increased proliferation in KO HSC/progenitors.
