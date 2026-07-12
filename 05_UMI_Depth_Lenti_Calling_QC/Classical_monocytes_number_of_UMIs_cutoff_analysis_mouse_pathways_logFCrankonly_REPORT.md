# Classical_monocytes_number_of_UMIs_cutoff_analysis_mouse_pathways_logFCrankonly.R

**Script:** `Classical_monocytes_number_of_UMIs_cutoff_analysis_mouse_pathways_logFCrankonly.R`

## What it does
Same as mouse_pathways.R but changes fgsea ranking to raw logFC (no -log10P weighting) and raises N_PERM to 5000, with a duplicated set.seed(42) call for reproducibility; retains cell-count plots, combined dendrogram heatmap, and master ComplexHeatmap sections.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, Lenti_detection_thresholds_SCTv2.csv.

## Outputs
GSEA_Monocytes_<comparison>_all_results.csv, GSEA_Monocytes_cell_counts_per_case.pdf/.png, GSEA_Monocytes_sensitivity_heatmap.pdf/.png, GSEA_Monocytes_combined_heatmap_FDR<x>.pdf/.png, GSEA_Master_Sensitivity_Mono_FDR0_05.pdf.

## Conceptual aim / target
Tests whether monocyte Hallmark pathway calls change when ranking genes by logFC alone instead of significance-weighted logFC, using the mouse-labeled msigdbr call.

## Note
Near-duplicate of mouse_pathways.R (ranking metric/N_PERM differ) and of human_pathways_logFCrankonly.R (only the msigdbr call syntax differs); the 6th member of a tightly related family of 6 files.
