# 10x WT vs Jak2V617F: compartment-matched proliferation comparison

**Script:** `c_join.py`

## What it does
Joins pre-extracted per-cell gene panels for WT and Jak2V617F 10x samples, regresses out cell-cycle score (for clustering purposes only), runs PCA + KMeans clustering, labels each cluster by its dominant lineage signature, and defines an "LT-HSC+MPP1" compartment as the HSC-signature-dominant clusters. Statistically compares %cycling (chi-square) and proliferation score (Mann-Whitney) between WT and Jak2V617F within that compartment, within "committed" cells, and pooled across all cells; also runs a paired pseudobulk t-test and draws a 3-panel PCA figure (clusters / genotype / compartment).

## Inputs
tenx_panel.csv, tenx_WT_panel.npy/tenx_Jak2V617F_panel.npy and their _obs.csv files — all produced upstream by t_panel.py.

## Outputs
tenx_joint_obs.csv, tenx_cluster_signatures.csv, tenx_LTHSC_MPP1_comparison.csv, fig5_tenx_clusters.png.

## Conceptual aim / target
Answer the project's central question for the Jak2V617F genotype: is proliferation specifically increased within the primitive LT-HSC+MPP1 compartment (not just in committed/differentiated cells), matching compartments rather than pooling the whole hierarchy.

## Conclusions / findings
This is the script whose LT-HSC+MPP1 result (WT vs Jak2V617F, chi2 p=5.5e-19, cycling up in mutant) is hardcoded into c_final.py's MASTER_matched_compartment.csv as the project's headline Jak2V617F finding.
