# combined_regions_analysis.R

**Script:** `combined_regions_analysis.R`

## What it does
Combines separately-sequenced "C-region" and "T-region" lentiviral clonal barcode tables into one unified wide table, then for each of four destination tissues computes per-mouse Euclidean and Bray-Curtis distances from progenitor-population/bone-marrow-site source samples to that destination, runs Friedman tests across sources, and plots ranked/faceted distance-to-destination plots.

## Inputs
C_regions_results/combined_barcodes_C_regions.tsv and T_region_result/combined_barcodes_T_region.tsv.

## Outputs
combined_barcodes_wide_CT.tsv; CT_<destination>_distances_per_mouse.csv, CT_<destination>_distance_friedman.csv, ranked/faceted plots per destination x metric.

## Conceptual aim / target
Determine which hematopoietic progenitor population and bone-marrow site is clonally closest to specific downstream tissue destinations, using lentiviral barcode data.

## Note
Distinct analysis — a clonal barcode-tracking/lineage-distance study, unrelated in method to the Lenti-ratio UMI-based classification scripts.
