# cMono_Seurat_Harmony_Lenti.R

**Script:** `cMono_Seurat_Harmony_Lenti.R`

## What it does
Same Seurat/Harmony pipeline as cMono_Seurat_Harmony.R but adds Lenti+/Lenti- classification (ratio >=4e-5 / <=9e-7) before downsampling, stores it as metadata, and generates Lenti-status UMAP plots plus HB-vs-LB and Lenti+-vs-Lenti- density-difference maps.

## Inputs
counts_mono (classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, expected pre-loaded).

## Outputs
mono_00_combined.pdf ... mono_12_lenti_per_sample.pdf, HB-vs-LB and Lenti+-vs-Lenti- density plots, classical_mono_seurat_final.rds.

## Conceptual aim / target
Test whether Lenti+ vs Lenti- monocytes occupy distinct transcriptional subclusters/UMAP regions, in addition to HB vs LB.

## Note
Near-identical to the sibling cMono_Seurat_Harmony_Lenti.R in the parent folder — differs only by a missing initial data-load block and minor formatting.
