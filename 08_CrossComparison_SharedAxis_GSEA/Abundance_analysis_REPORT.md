# Abundance_analysis.R

**Script:** `Abundance_analysis.R`

## What it does
Loads per-cluster cell counts from a Seurat/scRNA clustering run (26 clusters, labeled by cell type e.g. cMono, LT-HSC, GMPs), converts to per-sample proportions, and excludes progenitor/non-myeloid clusters. Builds a reusable bar+point plot function with significance stars, then runs four comparisons: HB vs LB within LentiPos (paired), HB vs LB within LentiNeg (paired), LentiPos vs LentiNeg within HB (unpaired), and same within LB, each with BH-adjusted t-tests. A second section correlates total cells per cluster against lentiviral UMI counts per cluster.

## Inputs
number_of_cells_per_cluster_for_each_sampleLenti_res0.6.txt (cluster x sample cell-count matrix); hardcoded lenti_umis vector per cluster.

## Outputs
HB_cell_proportions_LentiPos_vs_LentiNeg.pdf/png, LB_cell_proportions_LentiPos_vs_LentiNeg.pdf/png, LentiPos_cell_proportions_HB_vs_LB.pdf/png, LentiNeg_cell_proportions_HB_vs_LB.pdf/png, lenti_umi_correlation.pdf/png.

## Conceptual aim / target
Determine whether clonal burden (HB/LB) or lentiviral-construct status shifts the cellular composition of hematopoietic clusters, and whether cluster size predicts lentiviral marking coverage.

## Note
Distinct analysis (cell-proportion/composition, not DE) — the only cluster-composition script in the set.
