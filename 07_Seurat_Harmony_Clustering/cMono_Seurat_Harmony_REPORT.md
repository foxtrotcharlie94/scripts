# cMono_Seurat_Harmony.R

**Script:** `cMono_Seurat_Harmony.R`

## What it does
Seurat/Harmony clustering pipeline for classical monocytes: downsamples all cells to 4149 UMIs, builds Seurat object with sample/exposure/batch metadata, runs NormalizeData/PCA, Harmony batch correction, clusters (res 0.5), UMAP, and generates cluster/sample/exposure UMAP plots plus HB-vs-LB kernel-density difference maps.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt (as counts_mono, loaded upstream).

## Outputs
mono_00_combined.pdf through mono_09_pct_LB_density.pdf, classical_mono_seurat_final.rds.

## Conceptual aim / target
Establish a batch-corrected classical monocyte subcluster map and assess whether HB vs LB mice occupy different transcriptional regions.

## Note
Near-duplicate/predecessor of the two cMono_Seurat_Harmony_Lenti.R files — same pipeline but lacks Lenti+/Lenti- classification/plots.
