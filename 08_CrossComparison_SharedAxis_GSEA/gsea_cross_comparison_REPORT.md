# gsea_cross_comparison.R

**Script:** `gsea_cross_comparison.R`

## What it does
Compares a reference GSEA contrast (LentiNeg_HBvsLB) against three other pseudobulk contrasts, computing NES scatter/correlations, leading-edge Jaccard overlap, meta-GSEA (testing whether the reference's significant pathway/gene sets are enriched in each other comparison's ranking), and a "fraction captured" metric (cos^2, Deming slope) to quantify shared signal. An optional Part B re-fits edgeR on random data splits to estimate split-half NES reliability and rescale captured fraction by a noise ceiling.

## Inputs
Saved GSEA/DE CSVs under pop_dir (per-comparison folders with GSEA_*.csv, DE_*.csv); optionally count_matrix.txt and sample_and_clusterNumber_for_each_cell_res0.5.txt.

## Outputs
NES_correlations.csv, metaGSEA_results.csv, geneMetaGSEA_results.csv, captured_fraction.csv, reliability.csv, captured_fraction_ceiling_corrected.csv, plus PNGs (NES_scatter, metaGSEA_enrichment, captured_summary, captured_ceiling_corrected).

## Conceptual aim / target
Quantify how much of the lentiviral (Pos vs Neg) transcriptomic signature overlaps with the clone-burden (HB vs LB) signature, and whether that overlap is real vs within-experiment noise.

## Note
Closely related to shared_axis_pathway.R/shared_axis_heatmap.R (same pop_dir concept, different statistical framing) and orchestrated by run_missing_cross_comparison.R.
