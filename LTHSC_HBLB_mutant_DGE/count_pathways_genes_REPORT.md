# count_pathways_genes.R

**Script:** `count_pathways_genes.R`

## What it does
Builds summary bar plots of DEG counts (FDR<0.1, up/down) comparing "excl. HB3/HB4" vs "all samples" analyses, and separately counts significant GSEA pathways (FDR<0.05) per comparison x database (stacked and faceted by direction), for LT-HSCs and Classical Monocytes.

## Inputs
DEG (_limma_voom_results.csv) and GSEA (_GSEA_<db>.csv) files from LT-HSCs (1)/<comparison>[/results] and Classical_Monocytes (1)/....

## Outputs
DEG_counts_FDR0.1.pdf/png, pathway_counts_stacked.pdf/png, pathway_counts_by_database.pdf/png.

## Conceptual aim / target
Quantify how many genes/pathways are differentially expressed/enriched per comparison, and sensitivity to excluding two animals.

## Note
Near-duplicate of count_pathways_genes_beta.R — beta version uses different FDR handling and computes pathway counts for both sample-inclusion variants (this one only for excl HB3/HB4).
