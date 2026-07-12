# distance_analyses_all_metrics.R

**Script:** `distance_analyses_all_metrics.R`

## What it does
Generalized distance-to-destination analysis: for 3 datasets (C-only, T-only, CT-combined) x 4 destinations x 4 distance metrics (Euclidean, Bray-Curtis, cosine, Jensen-Shannon), computes per-mouse distances from each source to the destination, runs Friedman tests, and produces ranked/faceted plots.

## Inputs
combined_barcodes_C_regions.tsv, combined_barcodes_T_region.tsv, combined_barcodes_wide_CT.tsv.

## Outputs
Per dataset x destination: distances_per_mouse.csv, friedman.csv, ranked/faceted plots (PNG/PDF).

## Conceptual aim / target
Comprehensively test which progenitor source is clonally closest to each destination, and whether the ranking is robust across metric choice.

## Note
Generalizes combine_and_analyze_CT.R's STEP 2 to 4 metrics/3 datasets instead of 2 metrics on the CT-combined set only.
