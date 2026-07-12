# cluster_markers.R

**Script:** `cluster_markers.R`

## What it does
For five mutant-genotype Seurat objects (Tet2, Asxl1, Jak2, Dnmt3a, Dnmt3a_R878H) and two clustering resolutions, runs FindAllMarkers (with per-cluster fallback), FindConservedMarkers across genotype, and a canonical-gene DotPlot to help assign lineage identities.

## Inputs
<name>_seurat.rds objects under seurat_harmony_output/<dataset>/, produced by run_seurat_harmony_all_mutants.R.

## Outputs
<name>_markers_res<r>_all.csv, _top15.csv, _conserved_markers_res<r>.csv, _canonical_dotplot_res<r>.png per dataset/resolution.

## Conceptual aim / target
Annotate the unsupervised HSPC clusters from the 5-genotype proliferation study with hematopoietic lineage identities.

## Note
Superseded/duplicated by fix_cluster_markers.R, which redoes this job for 3 datasets because the original FindAllMarkers call silently dropped large progenitor clusters; finish_asxl1_res1.2_markers.R patches the remaining Asxl1 gap.
