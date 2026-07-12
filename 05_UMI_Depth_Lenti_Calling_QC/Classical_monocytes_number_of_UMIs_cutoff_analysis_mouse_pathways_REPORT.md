# Classical_monocytes_number_of_UMIs_cutoff_analysis_mouse_pathways.R

**Script:** `Classical_monocytes_number_of_UMIs_cutoff_analysis_mouse_pathways.R`

## What it does
Same 8-case, 4-comparison monocyte GSEA sensitivity pipeline as the human_pathways.R file, using msigdbr(species="Mus musculus", collection="MH", db_species="MM") for Hallmark gene sets and ranking fgsea by sign(logFC)*-log10(PValue) with N_PERM=1000; includes cell-count bar plots, combined heatmap-with-dendrogram, and final master ComplexHeatmap.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, Lenti_detection_thresholds_SCTv2.csv.

## Outputs
GSEA_Monocytes_<comparison>_all_results.csv, GSEA_Monocytes_cell_counts_per_case.pdf/.png, GSEA_Monocytes_sensitivity_heatmap.pdf/.png, GSEA_Monocytes_combined_heatmap_FDR<x>.pdf/.png, GSEA_Master_Sensitivity_Mono_FDR0_05.pdf.

## Conceptual aim / target
Same UMI-filtering sensitivity question as the other monocyte GSEA scripts, using the (correctly labeled) mouse Hallmark gene-set API call.

## Note
Near-duplicate of human_pathways.R — functionally identical except the msigdbr call syntax used to fetch the same mouse Hallmark gene sets differs.
