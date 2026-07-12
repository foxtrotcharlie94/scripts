# all_genes_comparison.R

**Script:** `all_genes_comparison.R`

## What it does
Two-part GSEA workflow on human CD14+ monocyte DGE tables: (1) auto-detects DESeq2/MAST format DGE files, computes two ranking metrics, and runs fgsea against Hallmark and GO:BP for every file found; (2) a visualization block for four TET2 CHIP/CH contrasts producing dotplots, barplots, cross-contrast NES heatmaps, and running-enrichment plots.

## Inputs
DGE_tables/DGE_tables/ folder (recursively scanned DESeq2/MAST tables) for four TET2 CHIP/CH vs Control comparisons.

## Outputs
Per-file GSEA CSVs, GSEA_all_results.csv, plus dotplot/barplot/heatmap/enrichment figures.

## Conceptual aim / target
Characterize Hallmark/GO:BP pathway enrichment in human CD14+ monocytes between TET2-mutant/CHIP/CH states and controls.

## Note
Distinct analysis (human bulk GSEA on external clinical cohort) — same underlying logic/content as DGE_tables/run_gsea.R + plot_gsea_human.R combined into one file.
