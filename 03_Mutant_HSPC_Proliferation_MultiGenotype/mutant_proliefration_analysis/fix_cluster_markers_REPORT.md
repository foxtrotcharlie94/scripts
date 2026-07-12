# fix_cluster_markers.R

**Script:** `fix_cluster_markers.R`

## What it does
Regenerates cluster markers for Tet2, Dnmt3a_R878H, and Asxl1 because earlier marker files were truncated/broken (default logfc.threshold=0.25 caused large progenitor clusters to return zero marker rows). Uses a low threshold + per-cluster fallback, downsamples cells per cluster, restricts to variable genes for speed, and uses atomic write-then-verify-then-rename CSV writes to avoid truncated output.

## Inputs
<name>_seurat.rds files in seurat_harmony_output/ for Tet2, Dnmt3a_R878H, Asxl1.

## Outputs
Corrected <name>_markers_res<r>_all.csv / _top15.csv (optionally _conserved_markers_res<r>.csv), plus canonical dotplot PNGs.

## Conceptual aim / target
Fix a systematic bug that made earlier per-cluster marker tables incomplete, ensuring every cluster gets ranked marker genes.

## Note
Near-duplicate/bugfix version of cluster_markers.R for a 3-dataset subset, with added speed/robustness engineering.
