# concordance_analysis.R

**Script:** `concordance_analysis.R`

## What it does
Compares DEG and GSEA results between "all samples" and "excluding HB3/HB4" analysis variants for four Lenti/HB-LB comparisons, for both LT-HSCs and Classical Monocytes. Produces gene-level logFC scatter plots (colored by significance overlap) and pathway-level GSEA NES scatter plots, annotating correlations and labeling top discordant genes/pathways.

## Inputs
_limma_*results.csv DEG tables and _GSEA_<db>.csv pathway tables from allSamples_analysis/ and excluding_H3_H4/ subfolders.

## Outputs
Combined multi-panel concordance_logFC.pdf/png and concordance_NES_<db>.pdf/png per dataset.

## Conceptual aim / target
Assess how robust the DEG/pathway findings are to inclusion/exclusion of two possibly outlier/confounding animals.

## Note
Distinct analysis, shares comparison/database definitions with count_pathways_genes* and gene_and_pathway_analysis_beta* scripts.
