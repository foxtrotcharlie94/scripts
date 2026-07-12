# cytotrace2_LT-HSCs.R

**Script:** `cytotrace2_LT-HSCs.R`

## What it does
Runs CytoTRACE2 on a Seurat object of LT-HSCs to score differentiation state, plots UMAPs/violins by HB/LB and LentiPos/Neg group; re-runs pooled Slingshot pseudotime and builds a combined CytoTRACE2+Slingshot figure; computes multiple effect-size framings (raw diff, %range, Cohen's d) for HB vs LB and LentiPos vs LentiNeg; builds publication-ready multi-panel figures.

## Inputs
LT_HSC_sub5000_seurat.rds from Sub5000_output/.

## Outputs
CytoTRACE_UMAP_score/relative.png, CytoTRACE_violin_group/sample.png, CytoTRACE_vs_pseudotime.png, CytoTRACE_stats.txt, updated seurat.rds, Combined_figure.png, CytoTRACE_effect_size_summary.csv, Fig1/Fig2 figures.

## Conceptual aim / target
Determine whether high clonal burden and/or Lenti-marking status shifts LT-HSCs toward a more differentiated transcriptional state.

## Note
Distinct analysis using the third-party CytoTRACE2 R package.
