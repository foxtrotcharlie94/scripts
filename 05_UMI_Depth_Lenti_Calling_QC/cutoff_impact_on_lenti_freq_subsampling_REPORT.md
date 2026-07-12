# cutoff_impact_on_lenti_freq_subsampling.R

**Script:** `cutoff_impact_on_lenti_freq_subsampling.R`

## What it does
Loads per-cell gene UMI matrices for LT-HSCs and Classical Monocytes, computes the original %Lenti+ per sample at a fixed cutoff (1e-4), then repeats after downsampling each cell's UMIs to several fixed depths (1000-15000), excluding samples below a minimum cell count; produces heatmaps and line/trend plots by exposure group.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt, classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt.

## Outputs
lenti_1e4_heatmap_hsc/mono.pdf, lenti_1e4_line_combined.pdf, lenti_1e4_persample_hsc/mono.pdf, lenti_1e4_subsampling_results.csv.

## Conceptual aim / target
Test whether the fraction of cells called Lenti-positive is a sequencing-depth artifact, by subsampling UMIs to a common depth.

## Note
Distinct analysis — single-cell UMI-depth QC/robustness check.
