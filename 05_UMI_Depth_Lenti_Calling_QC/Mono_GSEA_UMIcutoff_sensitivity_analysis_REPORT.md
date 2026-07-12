# Mono_GSEA_UMIcutoff_sensitivity_analysis.R

**Script:** `Mono_GSEA_UMIcutoff_sensitivity_analysis.R`

## What it does
GSEA sensitivity analysis for classical monocytes testing whether Hallmark pathway results depend on UMI-filtering strategy: 4 DE comparisons x 4 UMI-filtering cases (all cells, UMI-matched, fixed threshold >=13700, fixed threshold + matched), pseudobulk edgeR DE then fgsea against Hallmark, assembled into a 4x2 NES heatmap grid at two FDR cutoffs.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt and Lenti_detection_thresholds_SCTv2.csv.

## Outputs
GSEA_Monocytes_<comparison>_all_results.csv per comparison, combined GSEA_Monocytes_sensitivity_heatmap.pdf/.png.

## Conceptual aim / target
Test whether Hallmark pathway enrichment results in classical monocytes are robust to UMI-count filtering/matching choices.

## Note
Distinct but thematically related to the Classical_monocytes_number_of_UMIs_cutoff_analysis* family (simpler, 4-case version).
