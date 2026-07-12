# LTHSC_Seurat_Harmony.R

**Script:** `LTHSC_Seurat_Harmony.R`

## What it does
Loads the LT-HSC gene x cell UMI count matrix, downsamples every cell to 5000 UMIs, builds a Seurat object with sample/exposure/batch metadata, runs NormalizeData/PCA, Harmony batch-integration, clusters (res 0.5) and UMAP. Generates cluster/sample/exposure UMAP plots, per-sample highlight projections, and 2D kernel-density HB-vs-LB difference maps (raw diff, log2FC capped, %LB composition). Saves the Seurat object as RDS.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt.

## Outputs
00_combined.pdf through 09_pct_LB_density.pdf, LT-HSC_seurat_final.rds.

## Conceptual aim / target
Cluster LT-HSCs by transcriptome (batch-corrected) and visualize whether HB vs LB clones occupy distinct UMAP regions.

## Note
Precursor/simpler version of LTHSC_Seurat_Harmony_Lenti.R — same pipeline for LT-HSCs but without Lenti+/Lenti- classification/plots.
