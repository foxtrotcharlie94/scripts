# No_of_LentiAll_UMI_analysis.R

**Script:** `No_of_LentiAll_UMI_analysis.R`

## What it does
QC analysis of the LentiPos/LentiNeg ratio-threshold cutoff per cluster: scatter of log-total vs log-Lenti UMIs with cutoff lines, sequencing-depth-by-genotype violin plots, Lenti-status-fraction-vs-depth-bin curves, reference curves for what raw UMI count each ratio threshold corresponds to, distribution of LentiAll UMIs in LentiPos cells, hypergeometric-downsampling robustness of LentiPos calls, and cell counts surviving various UMI cutoffs.

## Inputs
lenti_ratio_per_cluster/per_cell_lenti_ratio.csv.

## Outputs
per_cluster_lentiStatus_QC_summary.csv and eight numbered PDF plots, plus downsampling-robustness CSVs, in lenti_cutoff_QC/.

## Conceptual aim / target
Verify the ratio-based cutoff is not a sequencing-depth artifact, and quantify robustness of individual cell calls to UMI subsampling.

## Note
Companion to abundance_analysis_ratioThr.R (same input/thresholds) — that script does the biological abundance analysis, this is the depth-bias/robustness QC. Byte-identical copy also exists at Downloads root as No_of_LentiAll_UMI_analysis.R.
