# DE_GSEA_analysis.R

**Script:** `DE_GSEA_analysis.R`

## What it does
A large (~700 line) cross-comparison diagnostics script comparing a reference contrast (LentiNeg_HBvsLB) against three others (LentiPos_HBvsLB, HB_PosVsNeg, LB_PosVsNeg) using saved GSEA CSVs: computes NES scatter plots with Pearson/Spearman correlation, leading-edge Jaccard overlap for jointly-significant pathways, "meta-GSEA" testing at multiple significance thresholds, and a "fraction captured" metric (cosine-similarity-based shared-signal / Deming-slope / residual) distinguishing "burden axis" vs "lenti axis" signal (Part A). It repeats an analogous gene-level meta-GSEA (Part A2), and has an optional Part B that re-derives split-half NES reliability from raw counts to rescale the captured-fraction into a noise-ceiling-corrected metric.

## Inputs
GSEA_<collection>_<comparison>_rankBy_<ranking>.csv and DE_<comparison>_byPValue.csv per comparison under pop_dir; Part B additionally reads count_matrix.txt and sample_and_clusterNumber_for_each_cell_res0.5.txt.

## Outputs
NES_correlations.csv, metaGSEA_results.csv, captured_fraction.csv, geneMetaGSEA_results.csv, leadingEdge_jaccard_*.csv, NES_scatter_*.png, metaGSEA_enrichment_*.png, geneMetaGSEA_enrichment_*.png, captured_summary.png, and (if Part B run) reliability.csv, captured_fraction_ceiling_corrected.csv, captured_ceiling_corrected.png.

## Conceptual aim / target
Determine how much of the "lentiviral-positive vs negative" transcriptional signature overlaps with the "high-burden vs low-burden" signature, both at the pathway and gene level, and how much of that overlap is real signal vs measurement noise.

## Note
Distinct, more elaborate analysis than Cross-contrast pathway heatmap.R — computes formal statistical overlap/reliability metrics rather than just a static union heatmap.
