# gene_and_pathway_analysis_beta_LTHSC.R

**Script:** `gene_and_pathway_analysis_beta_LTHSC.R`

## What it does
Runs limma-voom DE (5 pairwise HB/LB x LentiPos/LentiNeg comparisons) on LT-HSC count data excluding HB3/HB4, then GSEA (GO BP/MF/CC, KEGG, Hallmark) per comparison, filtering out non-hematopoietic tissue terms, producing volcano plots and combined GSEA dot plots.

## Inputs
DEG_* raw count/DEG tables per comparison subfolder under LT-HSCs (1)/.

## Outputs
Per comparison: _limma_voom_results.csv, _volcano.pdf/png, _GSEA_<db>.csv/plots, in a results/ subfolder.

## Conceptual aim / target
Identify genes/pathways differentially expressed between LentiPos/LentiNeg and HB/LB in LT-HSCs, excluding suspected outlier animals.

## Note
Near-duplicate of gene_and_pathway_analysis_beta_LT_HSCs_exludingH3H4.R (same exclusion, different DE method: voom/edgeR here vs limma-on-CPM there); parallels gene_and_pathway_analysis_beta_cMono.R (same code, cMono instead).
