# Classical_monocytes_number_of_UMIs_cutoff_analysis.R

**Script:** `Classical_monocytes_number_of_UMIs_cutoff_analysis.R`

## What it does
Runs a GSEA sensitivity analysis on Classical Monocytes across 4 group comparisons (LentiPos vs LentiNeg in LB, LentiPos vs LentiNeg in HB, HB vs LB in LentiNeg, HB vs LB in LentiPos) crossed with 6 nUMI-filtering "cases" (no filter, fixed threshold 13700, pooled NB threshold 9227, each with/without UMI-matching), each followed by pseudobulk edgeR DE and fgsea against MSigDB Hallmark (mouse). The script assembles per-comparison heatmaps at two FDR cutoffs (0.25/0.10) into a 4x2 grid, then appends a second block using ComplexHeatmap to build an "Ultra HD" master clustered heatmap with FDR-masked NES values.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt (gene x cell UMI matrix), Lenti_detection_thresholds_SCTv2.csv (per-sample thresholds); relies on a pre-existing `all_combined_res`/`col_fun`/`col_meta` object for the second block (not self-contained).

## Outputs
GSEA_Monocytes_<comparison>_all_results.csv per comparison, GSEA_Monocytes_sensitivity_heatmap.pdf/.png (4x2 grid), GSEA_Master_Sensitivity_Ultra_HD.pdf.

## Conceptual aim / target
Test how robust the Hallmark pathway calls for burden (HB/LB) and lentiviral-marking (LentiPos/Neg) contrasts in monocytes are to different UMI-depth filtering/matching choices.

## Note
Near-duplicate of Classical_monocytes_number_of_UMIs_cutoff_analysis_human_pathways.R and _mouse_pathways.R (same 6-case core loop); this file additionally appends the "Ultra HD" ComplexHeatmap block (which references undefined objects, i.e. is meant to run after another script) that the other two variants execute in fuller, self-contained form with 8 cases.
