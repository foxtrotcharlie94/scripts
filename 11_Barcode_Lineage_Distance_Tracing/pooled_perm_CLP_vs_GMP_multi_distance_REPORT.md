# pooled_perm_CLP_vs_GMP_multi_distance.R

**Script:** `pooled_perm_CLP_vs_GMP_multi_distance.R`

## What it does
Same DARLIN dataset and CLP-vs-GMP-to-T-lineage question as the bootstrap script, but runs a permutation test: pools CLP+GMP barcode mass, randomly reassigns per the observed global CLP fraction (5000 perms), computing a null distribution of delta for BH-adjusted p-values per bone x destination x metric.

## Inputs
combined_barcodes_C_regions.tsv and combined_barcodes_T_regions.tsv (same as bootstrap script).

## Outputs
pooled_perm_CLP_vs_GMP_multi_distance_results.csv (+wide version), tornado and raw-distance plots.

## Conceptual aim / target
Same as pooled_bootstrap_CLP_vs_GMP_destfixed.R — test whether T-lineage destinations are significantly closer to CLP or GMP — via permutation-based null instead.

## Note
Near-duplicate/companion of pooled_bootstrap_CLP_vs_GMP_destfixed.R — differs only in statistical inference method.
