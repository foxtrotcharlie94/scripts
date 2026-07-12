# gene_and_pathway_analysis_beta_cMono_exludingH3H4.R

**Script:** `gene_and_pathway_analysis_beta_cMono_exludingH3H4.R`

## What it does
Same limma-on-log-CPM DGE + GSEA pipeline as the cMono "_all_samples" CPM variant, but excludes HB3/HB4, writing to an excluding_H3_H4 output directory.

## Inputs
DEG_* CPM-style tables per comparison subfolder under Classical_Monocytes (1)/.

## Outputs
_limma_voom_results.csv, _volcano, _GSEA_<db>.csv/plots under Classical_Monocytes (1)/excluding_H3_H4/<comparison>/.

## Conceptual aim / target
HB/LB x LentiPos/Neg DEG/pathway analysis in classical monocytes, CPM-input limma method, outlier-excluded robustness arm.

## Note
Near-duplicate of gene_and_pathway_analysis_beta_cMono_all_samples.R; cell-type counterpart of gene_and_pathway_analysis_beta_LT_HSCs_exludingH3H4.R.
