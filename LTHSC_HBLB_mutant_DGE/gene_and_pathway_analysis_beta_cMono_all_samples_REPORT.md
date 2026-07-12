# gene_and_pathway_analysis_beta_cMono_all_samples.R

**Script:** `gene_and_pathway_analysis_beta_cMono_all_samples.R`

## What it does
Same limma-on-log-CPM DGE + GSEA pipeline as the LT-HSC CPM "_all_samples" variant, but for Classical Monocytes with ALL samples included, writing to an allSamples_analysis directory.

## Inputs
DEG_* CPM-style tables per comparison subfolder under Classical_Monocytes (1)/.

## Outputs
_limma_voom_results.csv, _volcano, _GSEA_<db>.csv/plots under Classical_Monocytes (1)/allSamples_analysis/<comparison>/.

## Conceptual aim / target
HB/LB x LentiPos/Neg DEG/pathway question for classical monocytes, all-samples arm, CPM-input limma method.

## Note
Cell-type counterpart of gene_and_pathway_analysis_beta_LT_HSCs_all_samples.R; parallels gene_and_pathway_analysis_beta_cMono_exludingH3H4.R.
