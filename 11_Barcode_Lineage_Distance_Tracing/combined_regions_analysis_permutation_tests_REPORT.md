# combined_regions_analysis_permutation_tests.R

**Script:** `combined_regions_analysis_permutation_tests.R`

## What it does
Self-contained pipeline that first combines separately-sequenced C-region and T-region barcode wide tables into one CT table, then pools UMIs across mice per sample_type and runs a UMI-level Binomial permutation test (5000 perms) comparing each hematopoietic source's Euclidean/Bray-Curtis distance to two destination cell populations (default Blood_T_cells vs Blood_Granulocytes).

## Inputs
combined_barcodes_C_regions.tsv and combined_barcodes_T_region.tsv from C_regions_results/ and T_region_result/.

## Outputs
combined_barcodes_wide_CT.tsv; per-source permutation result CSVs (Euclidean and Bray-Curtis); tornado bar plots showing significantly biased sources, BH-adjusted.

## Conceptual aim / target
Determine whether specific hematopoietic progenitor populations in specific bones show statistically significant clonal bias toward one blood destination over another.

## Note
Near-duplicate of files\pooled_perm_tests.R — an earlier/simpler version (2 distance metrics, single hardcoded contrast, inline CT-combine step) versus the generalized version (4 metrics, multiple contrasts, assumes CT file exists).
