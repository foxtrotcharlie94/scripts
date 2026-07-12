# H2B_analysis.R

**Script:** `H2B_analysis.R`

## What it does
Builds "label-retaining" (LR_up) and "label-losing" (LR_down) gene sets from the top/bottom 200 genes (by negated log2FoldChange) of an H2B-GFP-Low-vs-High HSC DEG reference, then runs clusterProfiler::GSEA of these two sets against 5 LT-HSC lentiviral/burden DEG comparisons. Produces per-comparison enrichment plots, a summary NES/FDR dotplot, leading-edge overlap heatmaps and Venn diagrams, ranked-list correlation scatterplots, plus a full H2B-ranked-list pathway analysis (GO:BP/KEGG/Hallmarks).

## Inputs
DEgene2_analysis_Res_H2BGFP_Low_VS_High.csv and 5 LT-HSC pseudobulk DEG text files under LT-HSCs (1)/<comparison>/.

## Outputs
GSEA_H2B_on_<comparison>.png (5 files), Summary_dotplot.png, LeadingEdge_heatmap/Venn PNGs, RankedList_scatter.png, plus GO:BP/KEGG/Hallmark bar/dot/csv under H2B_pathways/.

## Conceptual aim / target
Test whether the H2B label-retaining (quiescent) vs label-losing (dividing) HSC transcriptional signature is recapitulated in the lentiviral-marking and clonal-burden LT-HSC comparisons.

## Note
Near-duplicate of H2B_fdrcutoff\H2B_analysis.R — same pipeline/outputs; differs only in how LR_up/LR_down gene sets are defined (fixed top/bottom 200 here vs padj<0.05 in the other file).
