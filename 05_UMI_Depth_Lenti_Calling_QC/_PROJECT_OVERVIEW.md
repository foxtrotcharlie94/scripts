# Project overview — UMI-Depth & Lenti+/- Calling QC

Rigor/QC scripts testing whether the lentiviral-marking (Lenti+/Lenti-) cell classification and
downstream GSEA results are artifacts of sequencing depth (UMIs per cell) rather than true
biology: threshold-derivation methods (binomial, Beta-mixture, Wilson CI), UMI-matched/fixed-
threshold GSEA sensitivity analyses across many gene-set/ranking-metric variants for classical
monocytes and LT-HSCs, an inflammation-gene-panel heatmap repeated at multiple filtering
stringencies, and Lenti+ frequency-by-cluster analyses. Many of these scripts are large families
of near-duplicates differing only in a single filtering/ranking parameter (documented per-file in
each report) — this reflects genuine sensitivity-analysis practice (deliberately re-running the
same test many ways) rather than accidental duplication.
