# Classical_monocytes_number_of_UMIs_cutoff_analysis_human_pathways_logFCrankonly.R

**Script:** `Classical_monocytes_number_of_UMIs_cutoff_analysis_human_pathways_logFCrankonly.R`

## What it does
Identical structure to the human_pathways.R file (4 comparisons x 8 UMI cases, pseudobulk edgeR, cell-count plots, combined dendrogram heatmap, master ComplexHeatmap), but changes the fgsea ranking to raw logFC only (dropping the sign*-log10P weighting) and raises N_PERM to 5000; also sets set.seed(42) (called twice) before each fgsea call for reproducibility.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, Lenti_detection_thresholds_SCTv2.csv.

## Outputs
Same file set as human_pathways.R (GSEA_Monocytes_*_all_results.csv, cell-count plots, sensitivity heatmap, combined dendrogram heatmaps, GSEA_Master_Sensitivity_Mono_FDR0_05.pdf).

## Conceptual aim / target
Re-test the same UMI-filtering sensitivity question but using logFC-only gene ranking (rather than signed significance) to check whether pathway calls are sensitive to the ranking-metric choice.

## Note
Near-duplicate of human_pathways.R — only the fgsea ranking statistic and permutation count differ; also near-identical to the mouse_pathways_logFCrankonly.R version.
