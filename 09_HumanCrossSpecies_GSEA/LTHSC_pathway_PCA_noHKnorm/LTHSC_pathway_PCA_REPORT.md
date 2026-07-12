# LTHSC_pathway_PCA.R

**Script:** `LTHSC_pathway_PCA.R`

## What it does
Cross-dataset GSEA pipeline comparing six published Jakobsen human HSPC differential-expression tables (TET2-mutant vs WT, CH-WT vs non-CH, three aging cohorts) against five in-house LT-HSC pseudobulk DEG comparisons. Maps human genes to mouse orthologs via babelgene, builds up/down gene sets per dataset, runs GSEA against each LT-HSC ranked list, and additionally runs GO:BP/KEGG/Hallmark GSEA on the original human data.

## Inputs
Six Jakobsen CSV DE tables; five LT-HSC pseudobulk DEG .txt files under "LT-HSCs (1)/<comparison>/" folders.

## Outputs
Per-dataset folders under Jakobsen_HSPC_multi_GSEA/ with GSEA plots, Summary_dotplot.png, LeadingEdge heatmaps/Venn diagrams, RankedList_scatter.png, orthologs.csv, gene_sets.csv, plus global Combined_dotplot.png and Combined_summary.csv; pathway subfolder with GO/KEGG/Hallmark plots and CSVs.

## Conceptual aim / target
Determine whether human CH/TET2/aging HSPC transcriptional signatures are enriched in the mouse LT-HSC lentiviral-clone comparisons — cross-species validation of clonal hematopoiesis expression signatures.

## Note
Distinct analysis; no close duplicate among the R files.
