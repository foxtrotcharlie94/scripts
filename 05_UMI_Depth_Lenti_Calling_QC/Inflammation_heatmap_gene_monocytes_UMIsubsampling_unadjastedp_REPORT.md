# Inflammation_heatmap_gene_monocytes_UMIsubsampling_unadjastedp.R

**Script:** `Inflammation_heatmap_gene_monocytes_UMIsubsampling_unadjastedp.R`

## What it does
Essentially identical pipeline to Inflammation_heatmap_gene_monocytes_UMIsubsampling.R, but the significance mask on the logFC heatmap uses raw (unadjusted) PValue < 0.05 rather than FDR < 0.05.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, Lenti_detection_thresholds_SCTv2.csv.

## Outputs
inflammation_genes_DEG_equalUMI_pval_all_results.csv and Inflammation_genes_heatmap_Monocytes_equalUMI_pval.pdf.

## Conceptual aim / target
Same as the UMI-equalized inflammation analysis but with a more permissive significance threshold, as a sensitivity check.

## Note
Near-identical duplicate of Inflammation_heatmap_gene_monocytes_UMIsubsampling.R, differing only in FDR vs unadjusted-p mask.
