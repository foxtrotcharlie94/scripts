# gene_and_pathway_analysis_beta_LTHSC_all_samples.R

**Script:** `gene_and_pathway_analysis_beta_LTHSC_all_samples.R`

## What it does
Same limma-voom DGE + GSEA pipeline as beta_LTHSC.R, but includes ALL samples (no HB3/HB4 exclusion), writing to an allSamples_analysis output directory.

## Inputs
Same DEG_* raw count/DEG tables per comparison subfolder under LT-HSCs (1)/.

## Outputs
Same file types, written under LT-HSCs (1)/allSamples_analysis/<comparison>/.

## Conceptual aim / target
Same HB/LB x LentiPos/Neg DEG/pathway question as beta_LTHSC.R but as the "all samples" arm used for concordance/robustness checks.

## Note
Near-duplicate of gene_and_pathway_analysis_beta_LTHSC.R — identical except sample filter and output path.
