# Asxl1 10x: joined compartment-matched WT vs KO comparison (wk04 + wk36)

**Script:** `ax_join.py`

## What it does
Joins the 4 per-sample panel extracts written by ax_panel.py, cell-cycle-regresses (for clustering only) and PCA/MiniBatchKMeans clusters all cells together, labels clusters by lineage signature, defines the LT-HSC+MPP1 compartment, and compares WT vs KO %cycling (chi-square) and proliferation score (Mann-Whitney) within that compartment at each timepoint (wk04, wk36) and pooled, within whole LSK for comparison, plus a paired (across-timepoint) pseudobulk t-test.

## Inputs
/tmp/asxl1/panel.csv and the 4 sample *_X.npy/_obs.csv files written by ax_panel.py.

## Outputs
asxl1_10x_joint_obs.csv, asxl1_10x_comparison.csv, asxl1_10x_pseudobulk.csv; console cluster labels, compartment composition, and per-timepoint/pooled stats.

## Conceptual aim / target
Test whether Asxl1-KO shows increased proliferation specifically within the primitive LT-HSC+MPP1 compartment at two ages/timepoints (4 and 36 weeks), matching the same compartment-matched-comparison logic used for the other genotypes in this project.

## Conclusions / findings
Produces the batch-free (WT+KO co-captured), compartment-matched Asxl1 WT-vs-KO cycling/proliferation comparison at wk04, wk36 and pooled — this two-file (panel-extraction + join) design was later reused as the template for the Tet2, Dnmt3a and Jak2 per-genotype panel/join scripts below.
