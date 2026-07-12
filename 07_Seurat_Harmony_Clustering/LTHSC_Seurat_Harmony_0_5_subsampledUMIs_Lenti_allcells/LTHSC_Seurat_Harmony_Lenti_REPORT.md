# LTHSC_Seurat_Harmony_Lenti.R

**Script:** `LTHSC_Seurat_Harmony_Lenti.R`

## What it does
Extends LTHSC_Seurat_Harmony.R by additionally classifying each cell as Lenti+/Lenti-/Unclassified (ratio thresholds >=4e-5 / <=9e-7) before downsampling, then runs the same Seurat/Harmony/UMAP pipeline with lenti_status metadata added. Produces standard cluster/sample/exposure UMAP plots plus Lenti-status UMAP plots and generalized density-difference maps for both HB-vs-LB and Lenti+-vs-Lenti-.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt (including the "Lenti" pseudo-gene row).

## Outputs
hsc_00_combined.pdf through hsc_12_lenti_per_sample.pdf, plus HB-vs-LB and Lenti+-vs-Lenti- density plots.

## Conceptual aim / target
Determine whether Lenti+ and Lenti- cells occupy distinct transcriptional neighborhoods in LT-HSC UMAP space, in addition to the HB-vs-LB comparison.

## Note
Extended near-duplicate of LTHSC_Seurat_Harmony.R (adds Lenti+/- classification/plots); byte-identical copy also exists at the Downloads root (LTHSC_Seurat_Harmony_Lenti.R).
