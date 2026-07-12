# UMI_counts_per_Cell.R

**Script:** `UMI_counts_per_Cell.R`

## What it does
Computes HSP-family, Trp53, and mitochondrial UMI ratios per cell for LT-HSC and classical monocytes, then produces log-log scatter plots (Spearman correlation annotated) of HSP-vs-Trp53 and HSP-vs-Mitochondrial expression per cell, faceted by cell type and colored by exposure.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt, classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt.

## Outputs
hsp_vs_trp53_per_cell.pdf, hsp_vs_mt_per_cell.pdf, hsp_correlations_combined.pdf.

## Conceptual aim / target
Examine whether cellular stress markers co-vary at the single-cell level, to assess a shared stress/dissociation confound.

## Note
Companion to QC_Trp53_Mito_PCA_excl_Hsp.R — same three gene families, pairwise-correlation view rather than per-sample-distribution+PCA view.
