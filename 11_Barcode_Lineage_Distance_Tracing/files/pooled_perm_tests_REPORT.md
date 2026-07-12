# pooled_perm_tests.R

**Script:** `pooled_perm_tests.R`

## What it does
Reads the combined CT barcode table, pools UMIs across mice per (barcode, sample_type), and for several destination-pair contrasts runs a UMI-level Binomial permutation test (5000 perms) computing source-preference deltas under 4 distance metrics, producing tornado bar plots with BH-adjusted significance.

## Inputs
combined_barcodes_wide_CT.tsv.

## Outputs
CT_pooled_<pos>_vs_<neg>_perm_<metric>.csv and tornado plots per contrast x metric.

## Conceptual aim / target
Statistically test whether specific progenitor sources show significant clonal preference for one destination over another.

## Note
Generalization/superset of combined_regions_analysis_permutation_tests.R's permutation-test step (multiple contrasts, 4 metrics vs 2).
