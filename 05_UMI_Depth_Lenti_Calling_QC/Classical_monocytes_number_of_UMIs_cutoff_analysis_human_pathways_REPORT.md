# Classical_monocytes_number_of_UMIs_cutoff_analysis_human_pathways.R

**Script:** `Classical_monocytes_number_of_UMIs_cutoff_analysis_human_pathways.R`

## What it does
Same core pipeline as Classical_monocytes_number_of_UMIs_cutoff_analysis.R (4 comparisons x 8 UMI-filtering cases, including 2 new "sample-specific SCTv2 threshold" cases, pseudobulk edgeR + fgsea), but uses msigdbr(species = "Mus musculus", category = "H") gene sets and ranks fgsea by sign(logFC) * -log10(PValue). Adds a cell-count bar-plot section and a combined cross-comparison heatmap-with-dendrogram block, then a final ComplexHeatmap "Master Sensitivity" heatmap reading back its own written CSVs.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, Lenti_detection_thresholds_SCTv2.csv.

## Outputs
GSEA_Monocytes_<comparison>_all_results.csv, GSEA_Monocytes_cell_counts_per_case.pdf/.png, GSEA_Monocytes_sensitivity_heatmap.pdf/.png, GSEA_Monocytes_combined_heatmap_FDR<x>.pdf/.png, GSEA_Master_Sensitivity_Mono_FDR0_05.pdf.

## Conceptual aim / target
Same as the base file (UMI-filter sensitivity of Hallmark GSEA in monocytes), extended to 8 filtering cases and self-contained through to the master heatmap, despite filename saying "human_pathways" — gene sets used are actually mouse.

## Note
Near-duplicate of Classical_monocytes_number_of_UMIs_cutoff_analysis_mouse_pathways.R (identical except mouse_pathways version uses collection="MH", db_species="MM" msigdbr call instead of category="H"); also near-duplicate of the logFCrankonly variant but differs in fgsea ranking metric and N_PERM (1000 vs 5000).
