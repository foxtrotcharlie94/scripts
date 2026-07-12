# nUMIs_Lentipos_vs_Lentineg.R

**Script:** `nUMIs_Lentipos_vs_Lentineg.R`

## What it does
For LT-HSCs, runs a GSEA sensitivity analysis across 4 comparisons crossed with 4 UMI-filtering/matching cases and 2 FDR cutoffs, using pseudobulk edgeR + fgsea against Hallmark gene sets, assembled into a 4x2 grid of NES heatmaps.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt and Lenti_detection_thresholds_SCTv2.csv.

## Outputs
Per-comparison GSEA_LTHSCs_<id>_all_results.csv and combined GSEA_LTHSCs_sensitivity_heatmap.pdf/.png.

## Conceptual aim / target
Test whether GSEA/pathway results for HB-vs-LB and LentiPos-vs-LentiNeg in LT-HSCs are robust to nUMI-depth confounding.

## Note
UMI-depth robustness/sensitivity check, related in theme to gsea_cross_comparison.R but methodologically independent.
