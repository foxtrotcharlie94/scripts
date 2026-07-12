# pooled_bootstrap_CLP_vs_GMP_destfixed.R

**Script:** `pooled_bootstrap_CLP_vs_GMP_destfixed.R`

## What it does
Using DARLIN lentiviral-barcode data pooled across mice, computes Bray-Curtis/Cosine/Jensen-Shannon distances from two candidate progenitor sources (CLP, GMP) in each of 4 bones to two fixed destinations (Blood T cells, Thymus DN1), then bootstraps the source populations at observed sequencing depth to build a percentile CI on delta = distance(GMP,dest) - distance(CLP,dest).

## Inputs
combined_barcodes_C_regions.tsv and combined_barcodes_T_regions.tsv (DARLIN lineage-tracing dataset).

## Outputs
pooled_bootstrap_CLP_vs_GMP_destfixed_results.csv, tornado_CLP_vs_GMP_destfixed.png/pdf, raw_distance_CLP_vs_GMP_destfixed.png/pdf.

## Conceptual aim / target
Determine whether T-lineage cells derive preferentially from CLP or GMP progenitors, accounting for sampling noise.

## Note
Closely related to pooled_perm_CLP_vs_GMP_multi_distance.R — same question/data/metrics, source-side bootstrap CI here vs permutation-test null there.
