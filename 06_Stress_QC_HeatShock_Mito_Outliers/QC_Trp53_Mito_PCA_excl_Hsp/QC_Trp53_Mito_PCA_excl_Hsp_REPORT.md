# QC_Trp53_Mito_PCA_excl_Hsp.R

**Script:** `QC_Trp53_Mito_PCA_excl_Hsp.R`

## What it does
QC script for LT-HSCs examining per-cell Trp53, heat-shock-protein, and mitochondrial gene expression ratios via violin/boxplots per sample, then runs pseudobulk PCA on CPM-normalized sample aggregates twice — with all genes and excluding HSP+MT genes — to check whether stress-gene artifacts drive sample separation.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt.

## Outputs
hsc_qc_trp53/hsp/mt_per_sample.pdf, hsc_qc_combined.pdf, hsc_pseudobulk_pca_all_genes.pdf, hsc_pseudobulk_pca_no_hsp_mt.pdf, hsc_pseudobulk_pca_combined.pdf, hsc_qc_ratios_per_cell.csv.

## Conceptual aim / target
Assess whether stress/dissociation-artifact genes (Trp53/Hsp/mt) confound the pseudobulk PCA separation of HB vs LB samples in LT-HSCs.

## Note
Related to UMI_counts_per_Cell.R (same 3 gene families) but distinct — QC+PCA focused rather than pairwise correlation focused.
