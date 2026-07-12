# Lenti_frequency_in_subclusters_nalaysis.R

**Script:** `Lenti_frequency_in_subclusters_nalaysis.R`

## What it does
Computes Lenti+ cell frequency per sample and per Seurat cluster in LT-HSC and classical monocyte datasets: bar/stacked plots of Lenti+/Lenti-/Unclassified percentages, a Mono/HSC frequency ratio plot, per-cluster dotplots, a sanity check recomputing Lenti status from raw counts vs stored metadata, ratio-vs-absolute-UMI distribution plots, and a UMI-subsampling robustness analysis.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt, classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt, LT_HSC_seurat_final.rds, classical_mono_seurat_final.rds.

## Outputs
Many PDFs (lenti_freq_*_barplot/stacked, ratio_mono_over_hsc, cluster_dotplot_*, cutoff_comparison_*, distribution, subsampling_*_safe) plus matching CSVs.

## Conceptual aim / target
Characterize how Lenti+ clone frequency varies by sample, exposure, cell type, and subcluster, and test robustness to sequencing depth/cutoff definitions.

## Note
Near-duplicate of Lenti_frequency_in_subclusters_nalaysis_2.R — file _2 is a superset with two extra appended sections (HB-vs-LB statistical testing, UMI-depth distribution plots).
