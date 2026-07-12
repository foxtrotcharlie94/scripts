# finish_asxl1_res1.2_markers.R

**Script:** `finish_asxl1_res1.2_markers.R`

## What it does
Minimal targeted rerun that computes FindAllMarkers (per-cluster fallback, no conserved-marker step) for only the Asxl1 dataset at resolution 1.2, skipping the slow FindConservedMarkers step.

## Inputs
Asxl1_seurat.rds from seurat_harmony_output/Asxl1/.

## Outputs
Asxl1_markers_res1.2_all.csv and Asxl1_markers_res1.2_top15.csv.

## Conceptual aim / target
Quickly patch in the one missing piece (Asxl1 fine resolution) without re-running the entire slow pipeline.

## Note
A narrow, single-dataset/single-resolution variant of cluster_markers.R/fix_cluster_markers.R.
