# combine_barcodes.R

**Script:** `combine_barcodes.R`

## What it does
Generic utility that reads all per-sample two-column (barcode, UMI count) text files matching a pattern, derives sample names, sums duplicate rows, and pivots to a single wide barcode x sample matrix filled with 0 for absent barcodes.

## Inputs
Per-sample *_Barcode_UMI_number.txt files.

## Outputs
combined_barcodes_wide.tsv.

## Conceptual aim / target
Foundational data-wrangling step consolidating raw per-sample barcode-UMI files into the wide table used by all downstream barcode-lineage scripts.

## Note
Upstream utility producing inputs consumed by barcode_visualizations.R, distance_analyses_all_metrics.R, combine_and_analyze_CT.R, and the permutation-test scripts.
