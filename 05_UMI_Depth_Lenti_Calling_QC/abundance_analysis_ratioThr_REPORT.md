# abundance_analysis_ratioThr.R

**Script:** `abundance_analysis_ratioThr.R`

## What it does
Classifies cells as LentiPos/LentiNeg using ratio thresholds (>=5e-5 / <=1e-5), computes per-sample cluster-composition proportions across ~20 annotated Harmony clusters, and generates four significance-tested bar+dot plots (HB LentiPos vs LentiNeg, LB LentiPos vs LentiNeg, LentiPos HB vs LB, LentiNeg HB vs LB) plus a cluster-size vs Lenti-UMI correlation scatter.

## Inputs
lenti_ratio_per_cluster/per_cell_lenti_ratio.csv.

## Outputs
props_per_sample_genotype_cluster.csv, four PDF+PNG proportion plots, lenti_coverage_per_cluster.csv, lenti_umi_correlation.pdf/.png.

## Conceptual aim / target
Determine whether Lenti+ clones are over/under-represented in particular differentiation clusters, and whether this differs by burden.

## Note
Companion to No_of_LentiAll_UMI_analysis.R (same thresholds/input CSV) — this file focuses on abundance statistics, the other on UMI-depth QC for the same classification.
