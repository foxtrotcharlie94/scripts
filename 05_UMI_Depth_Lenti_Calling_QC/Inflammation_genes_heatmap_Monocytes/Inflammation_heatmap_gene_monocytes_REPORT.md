# Inflammation_heatmap_gene_monocytes.R

**Script:** `Inflammation_heatmap_gene_monocytes.R`

## What it does
Builds a curated ~100-gene inflammation panel and tests it across 4 comparisons x 8 UMI-filtering/matching cases in classical monocytes. For each combination pseudobulks per sample/group, runs edgeR QLF, and assembles a logFC matrix masked by FDR<0.05 into a ComplexHeatmap grouped by gene category and comparison.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt and Lenti_detection_thresholds_SCTv2.csv.

## Outputs
inflammation_genes_DEG_all_results.csv and Inflammation_genes_heatmap_Monocytes.pdf (logFC heatmap, FDR<0.05).

## Conceptual aim / target
Test the robustness of inflammatory gene expression differences (by lentiviral status and clone burden) in classical monocytes across UMI-depth filtering/matching strategies.

## Note
Near-duplicate of Inflammation_heatmap_gene_monocytes_UMIsubsampling.R — same gene panel/comparisons/cases/edgeR logic; this version has no UMI-equalization/subsampling step.
