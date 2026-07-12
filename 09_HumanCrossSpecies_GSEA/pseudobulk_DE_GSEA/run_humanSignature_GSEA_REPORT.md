# run_humanSignature_GSEA.R

**Script:** `run_humanSignature_GSEA.R`

## What it does
Builds human UP/DOWN gene signatures (top/bottom N by logFC among significant genes) from 4 human TET2/CHIP monocyte DGE files, maps mouse pseudobulk DE results (3 monocyte subsets x 4 comparisons) to human orthologs via babelgene, ranks genes two ways, and runs fgsea to test enrichment of each human signature in each mouse contrast; produces per-file enrichment PDFs and four cross-cutting summary figures.

## Inputs
4 human DGE files under DGE_tables/ and per-subset/comparison mouse DE_*.csv files under pseudobulk_DE_GSEA/<subset>/<comparison>/.

## Outputs
Per-run CSVs, humanSignature_GSEA_all_results.csv, and figures (enrichment PDFs, NES heatmap, UP-vs-DOWN concordance, cross-context correlation, metric agreement summaries).

## Conceptual aim / target
Test whether the mouse HB/LB and LentiPos/Neg monocyte transcriptional programs recapitulate known human TET2-mutant/CHIP monocyte expression signatures.

## Note
Cross-species translational extension of pseudobulk_DE_GSEA.R's mouse DE outputs; conceptually related to gsea_cross_comparison.R but uses human-derived reference signatures.
