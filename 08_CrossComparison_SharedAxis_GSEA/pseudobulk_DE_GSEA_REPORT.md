# pseudobulk_DE_GSEA.R

**Script:** `pseudobulk_DE_GSEA.R`

## What it does
Self-contained pseudobulk DE + GSEA pipeline: loads a genes x cells count matrix and cluster/sample metadata, classifies cells as LentiPos/LentiNeg via UMI ratio, pools clusters into a named population, and runs 4 edgeR pseudobulk comparisons each with PCA diagnostics, BCV/MD plots, volcano plots, sorted DE CSVs, and Hallmark + GO_BP fgsea (ranked by logFC and logFC x -log10P), finishing with group-level summary bar plots.

## Inputs
count_matrix.txt (genes x cells), sample_and_clusterNumber_for_each_cell_res0.5.txt, optionally barcode_counts_raw_with_sort.csv.

## Outputs
Per-comparison folders with PCA/BCV/MD/volcano PNGs, 3 sorted DE_*.csv variants, GSEA_*.csv/.png per collection x ranking; base-level summary_DE_counts.csv, summary_pathway_counts.csv, summary.png.

## Conceptual aim / target
Characterize how gene expression/pathway activity differ by lentiviral marking status and clonal burden within a chosen cell population — the core DE/GSEA engine of the whole clonal-hematopoiesis burden study.

## Note
The master/generic pipeline whose outputs (GSEA_*.csv, DE_*.csv) are the direct inputs to gsea_cross_comparison.R, shared_axis_heatmap.R, shared_axis_pathway.R.
