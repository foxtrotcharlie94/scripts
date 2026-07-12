# combine_and_analyze_CT.R

**Script:** `combine_and_analyze_CT.R`

## What it does
STEP 1 combines C-region and T-region barcode wide tables into one CT table. STEP 2 computes per-mouse Euclidean and Bray-Curtis distances from hematopoietic sources to four destinations, normalizing C and T barcodes separately then averaging 50/50; runs Friedman tests and produces ranked/faceted bar plots.

## Inputs
combined_barcodes_C_regions.tsv, combined_barcodes_T_region.tsv.

## Outputs
combined_barcodes_wide_CT.tsv; per-destination distances_per_mouse.csv and friedman.csv; ranked and faceted distance plots.

## Conceptual aim / target
Identify which bone-marrow progenitor population/location is clonally closest to each downstream destination tissue.

## Note
Near-duplicate/superset of distance_analyses_all_metrics.R's CT-dataset branch (that one generalizes to 4 metrics x 3 datasets).
