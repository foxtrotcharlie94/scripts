# Inflammation_heatmap_gene_monocytes_UMIsubsampling.R

**Script:** `Inflammation_heatmap_gene_monocytes_UMIsubsampling.R`

## What it does
Same inflammation-gene-panel x 4-comparison x 8-case pipeline as the sibling script, but adds a global-minimum-UMI-target subsampling step before pseudobulking so every group is built from an equalized total UMI budget before edgeR QLF is run.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, Lenti_detection_thresholds_SCTv2.csv.

## Outputs
inflammation_genes_DEG_equalUMI_all_results.csv and Inflammation_genes_heatmap_Monocytes_equalUMI.pdf.

## Conceptual aim / target
Same biological question as the non-subsampled version, additionally controlling for total sequencing depth per pseudobulk group.

## Note
Near-duplicate of Inflammation_heatmap_gene_monocytes.R (adds UMI-equalization) and of Inflammation_heatmap_gene_monocytes_UMIsubsampling_unadjastedp.R (identical except mask uses FDR vs raw p-value).
