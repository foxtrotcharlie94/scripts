# barcode_visualizations.R

**Script:** `barcode_visualizations.R`

## What it does
Reads a combined wide barcode/UMI table, parses sample names into mouse/population/source metadata, and generates a barcode-sharing histogram, per-sample/per-mouse barcode count bar charts, per-mouse UpSet plots, per-mouse-per-population Venn diagrams, and a heatmap of barcode counts by mouse x sample_type.

## Inputs
combined_barcodes_wide.tsv.

## Outputs
01_barcode_sharing_histogram, 02_barcodes_per_sample_per_mouse, 03a_upset_mouse<N>, 03b_venn_mouse<N>_<population>, 04_heatmap_all_mice PNGs/PDFs, plus barcodes_per_sample.csv and heatmap CSV.

## Conceptual aim / target
Descriptive/exploratory QC of lentiviral barcode-lineage tracing data.

## Note
Distinct analysis (descriptive visualization) feeding the same barcode-lineage project as combine_barcodes.R, distance_analyses_all_metrics.R, combine_and_analyze_CT.R, pooled_perm_tests.R.
