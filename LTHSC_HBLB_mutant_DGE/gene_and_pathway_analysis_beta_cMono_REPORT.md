# gene_and_pathway_analysis_beta_cMono.R

**Script:** `gene_and_pathway_analysis_beta_cMono.R`

## What it does
Same limma-voom DGE + GSEA pipeline as beta_LTHSC.R, but for Classical Monocytes, excluding HB3/HB4; produces volcano and GSEA dot plots per comparison.

## Inputs
DEG_* raw count/DEG tables per comparison subfolder under Classical_Monocytes (1)/.

## Outputs
<comp>_limma_voom_results.csv, _volcano.pdf/png, _GSEA_<db>.csv/plots under Classical_Monocytes (1)/<comparison>/results/.

## Conceptual aim / target
Identify DEGs/pathways between LentiPos/LentiNeg and HB/LB in classical monocytes, excluding possible outlier animals.

## Note
Cell-type counterpart of gene_and_pathway_analysis_beta_LTHSC.R — identical code, only base_dir changed.
