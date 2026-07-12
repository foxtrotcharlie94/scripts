# separate_regions_analysis.R

**Script:** `separate_regions_analysis.R`

## What it does
Multi-part concatenated script for DARLIN lineage-barcode data: combines per-sample barcode/UMI files into one wide matrix, produces barcode-sharing/overlap visualizations (histogram, bar plots, UpSet, Venn, heatmap), then four repeated destination-distance blocks (Thymus/Blood_T/Blood_Granulocytes/Heart_Granulocytes) each computing per-mouse Euclidean/Bray-Curtis distance from bone-marrow progenitor sources to a destination, with Friedman tests and ranked/faceted plots.

## Inputs
Per-sample *_Barcode_UMI_number.txt files under T_region_result/results_T_regions.

## Outputs
combined_barcodes_T_region.tsv; barcode-sharing/overlap PNGs/PDFs/CSVs; per-destination distances_per_mouse.csv, friedman.csv, ranked/faceted distance plots.

## Conceptual aim / target
Characterize clonal barcode-sharing patterns and statistically test which bone-marrow progenitor population/bone is the closest clonal source for each downstream tissue destination.

## Note
Shares the DARLIN dataset and CLP/GMP-vs-destination distance logic with the pooled_bootstrap/pooled_perm scripts, but is per-mouse and includes the barcode-combining/visualization steps; its four destination blocks are near-duplicates of each other.
