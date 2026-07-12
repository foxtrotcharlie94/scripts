# gene_and_pathway_analysis_beta_LT_HSCs_all_samples.R

**Script:** `gene_and_pathway_analysis_beta_LT_HSCs_all_samples.R`

## What it does
Same 5-comparison LT-HSC pipeline but uses limma directly on log2-CPM input (mean-log2CPM>0.5 filter) rather than voom/DGEList, with ALL samples included; runs the same 5-database GSEA and volcano/dot plots.

## Inputs
DEG_* CPM-style tables per comparison subfolder under LT-HSCs (1)/.

## Outputs
_limma_voom_results.csv (despite the CPM method), _volcano, _GSEA_<db>.csv/plots, under LT-HSCs (1)/allSamples_analysis/<comparison>/.

## Conceptual aim / target
Same HB/LB x LentiPos/Neg DEG/pathway question for LT-HSCs, all-samples arm, using a CPM-input limma model instead of voom.

## Note
Near-duplicate of gene_and_pathway_analysis_beta_LTHSC_all_samples.R (differs in DE method); parallel to _exludingH3H4 version (same method, HB3/HB4 excluded).
