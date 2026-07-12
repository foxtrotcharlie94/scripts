# cluser_composition_analysis.R

**Script:** `cluser_composition_analysis.R`

## What it does
Loads LT-HSC and classical monocyte Seurat objects, runs FindAllMarkers per cluster and saves top-50-marker tables plus top-10-marker DoHeatmap plots; a second section defines named cluster labels (e.g. "Quiescent LT-HSC", "Cycling HSC") and computes per-sample, per-cluster, per-(exposure x Lenti-status) group cell-frequency statistics with bar/dot/stacked-bar visualizations.

## Inputs
LT_HSC_seurat_final.rds, classical_mono_seurat_final.rds.

## Outputs
hsc/mono_cluster<N>_top50_markers.txt, cluster_markers_heatmap.pdf, combined_cluster_markers_heatmap.pdf, cluster_freq_barplot/dotplot/stacked.pdf, cluster_freq_per_sample.csv/_summary.csv.

## Conceptual aim / target
Identify marker genes defining each subcluster and quantify how the four groups (HB/LB x Lenti+/-) differ in distribution across annotated subclusters.

## Note
Distinct analysis — the marker-identification/cluster-composition-by-group step built on the Seurat objects from the cMono/LT-HSC Harmony scripts.
