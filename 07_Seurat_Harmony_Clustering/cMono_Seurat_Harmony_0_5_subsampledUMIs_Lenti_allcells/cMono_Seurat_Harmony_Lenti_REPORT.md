# cMono_Seurat_Harmony_Lenti.R

**Script:** `cMono_Seurat_Harmony_Lenti.R`

## What it does
Essentially identical to the "New folder" copy: classical monocyte Seurat/Harmony pipeline with Lenti+/Lenti- classification, Harmony batch correction, clustering/UMAP, Lenti-status UMAP plots, and HB-vs-LB plus Lenti+-vs-Lenti- density-difference maps — this version additionally loads counts_mono directly from file and prints a cluster-proportion table.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt (loaded directly at the top).

## Outputs
Same set as the sibling file.

## Conceptual aim / target
Same as sibling — classical monocyte subclustering plus Lenti+/- and HB/LB spatial comparison.

## Note
Near-duplicate of "New folder\cMono_Seurat_Harmony_Lenti.R" — same analysis, one being a slightly more complete/standalone version of the other.
