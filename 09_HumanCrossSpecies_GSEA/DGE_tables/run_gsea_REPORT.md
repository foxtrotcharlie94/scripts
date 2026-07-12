# run_gsea.R

**Script:** `run_gsea.R`

## What it does
Auto-discovers all DESeq2 (.tsv) and MAST (.csv) DGE result files recursively under INPUT_DIR, parses each into a common gene/logFC/padj format, builds two ranking metrics (raw logFC, and logFC x -log10(padj) with a zero-padj floor), and runs fgsea against human Hallmark and GO:BP gene sets for every file x metric x collection combination.

## Inputs
DESeq2_*/*.tsv and MAST_*/*.csv DGE tables under DGE_tables/DGE_tables/.

## Outputs
One <tag>__<metric>__<collection>.csv per file/metric/collection, plus a combined GSEA_all_results.csv, all written to DGE_tables/GSEA_results/.

## Conceptual aim / target
Systematically compute GSEA pathway enrichment for every available human monocyte DGE comparison (TET2 CHIP/CH/mutant contrasts) so downstream plotting scripts can visualize/compare them.

## Note
Distinct analysis — this is the upstream GSEA-computation engine that plot_gsea_human.R consumes.
