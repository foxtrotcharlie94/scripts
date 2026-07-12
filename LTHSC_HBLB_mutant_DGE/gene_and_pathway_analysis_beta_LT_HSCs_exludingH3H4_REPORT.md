# gene_and_pathway_analysis_beta_LT_HSCs_exludingH3H4.R

**Script:** `gene_and_pathway_analysis_beta_LT_HSCs_exludingH3H4.R`

## What it does
Same limma-on-log-CPM DGE + GSEA pipeline as the "_all_samples" CPM variant, but explicitly excludes HB3/HB4, writing to an excluding_H3_H4 results directory.

## Inputs
DEG_* CPM-style tables per comparison subfolder under LT-HSCs (1)/.

## Outputs
_limma_voom_results.csv, _volcano, _GSEA_<db>.csv/plots, under LT-HSCs (1)/excluding_H3_H4/<comparison>/.

## Conceptual aim / target
HB/LB x LentiPos/Neg DEG/pathway analysis in LT-HSCs using the CPM-input limma method, restricted to the outlier-excluded robustness arm.

## Note
Near-duplicate of gene_and_pathway_analysis_beta_LT_HSCs_all_samples.R — differs only in the HB3/HB4 exclusion filter and output directory.
