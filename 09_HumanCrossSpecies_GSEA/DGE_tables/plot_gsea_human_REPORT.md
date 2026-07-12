# plot_gsea_human.R

**Script:** `plot_gsea_human.R`

## What it does
Builds visualization outputs (dotplots, barplots, cross-contrast NES heatmaps, running-ES enrichment plots) for 4 specific human CD14+ monocyte TET2 CHIP/CH DGE comparisons, reading the GSEA result CSVs produced by run_gsea.R for the first three plot types, and recomputing ranked lists + re-fetching Hallmark/GO:BP gene sets on the fly to draw enrichment curves for the top significant pathways.

## Inputs
GSEA_results/*.csv (from run_gsea.R) plus the four raw DGE files (DESeq2 CHIP-vs-Control and CH-vs-Control TSVs, MAST Mutant-vs-Control/Wildtype CSVs).

## Outputs
dotplot_<collection>.pdf, barplot_<collection>.pdf, heatmap_<collection>_<metric>.pdf, enrichment_<label>.pdf, under GSEA_results/figures/.

## Conceptual aim / target
Visualize which Hallmark/GO:BP pathways are enriched/depleted in human TET2-mutant/CHIP monocytes relative to controls, across two DE methods (DESeq2, MAST) and two ranking metrics.

## Note
Downstream consumer of run_gsea.R output; distinct analysis (human data) from the mouse-focused monocyte/HSC scripts.
