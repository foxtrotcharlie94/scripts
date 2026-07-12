# count_pathways_genes_beta.R

**Script:** `count_pathways_genes_beta.R`

## What it does
Updated version of count_pathways_genes.R: computes DEG counts comparing all-samples vs excl-HB3/HB4, and pathway counts (FDR<0.05) for BOTH analysis variants, stacked by database and faceted by up/down direction.

## Inputs
_limma_*results.csv and _GSEA_<db>.csv from allSamples_analysis/ and excluding_H3_H4/ subfolders.

## Outputs
DEG_counts_FDR0.1.pdf/png and per-variant pathway_counts_stacked/_by_database.pdf/png.

## Conceptual aim / target
Same as count_pathways_genes.R but more complete — summarize DEG/pathway burden across both inclusion strategies.

## Note
Refined near-duplicate of count_pathways_genes.R.
