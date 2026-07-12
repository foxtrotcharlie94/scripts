# run_seurat_harmony_all_mutants.R

**Script:** `run_seurat_harmony_all_mutants.R`

## What it does
End-to-end Seurat pipeline for five public mouse HSPC scRNA-seq datasets (Tet2 GSE209994, Asxl1 GSE158184, Jak2 GSE227026, Dnmt3a GSE272266, Dnmt3a-R878H GSE233963): loads/demultiplexes raw 10x matrices (incl. HTO demux), QC filtering, cell-cycle scoring/regression, PCA, Harmony batch correction, UMAP, clustering at two resolutions, then marker/conserved-marker finding and canonical DotPlots, saving a full Seurat .rds per dataset.

## Inputs
Raw 10x mtx/h5/tar files per GEO accession under the mutant_proliefration_analysis base folder.

## Outputs
Per dataset: UMAP PNGs, metadata CSV, marker/conserved-marker CSVs, canonical dotplots, cluster x genotype composition CSVs, saved <name>_seurat.rds.

## Conceptual aim / target
Build a harmonized, clustered multi-genotype HSPC atlas to compare proliferation/lineage-composition phenotypes across CH driver mutations.

## Note
Primary/original pipeline producing the .rds objects consumed by cluster_markers.R, fix_cluster_markers.R, finish_asxl1_res1.2_markers.R.
