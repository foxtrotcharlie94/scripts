# Smart-seq2 pipeline step 2: Leiden clustering + per-cluster proliferation

**Script:** `c2_cluster_ss2.py`

**Status:** Pipeline step 2 of 2 (paired with c1_build_ss2.py)

## What it does
Loads the AnnData built by c1_build_ss2.py, computes a neighbor graph, Leiden clusters and UMAP embedding, labels each cluster by its dominant lineage-marker signature, and computes %cycling and mean proliferation score per (cell-type cluster x genotype) combination.

## Inputs
ss2_built.h5ad (output of c1_build_ss2.py).

## Outputs
ss2_cluster_proliferation.csv, ss2_cluster_marker_scores.csv, fig4_ss2_umap.png, ss2_clustered.h5ad.

## Conceptual aim / target
Determine whether Tet2/Dnmt3a-associated proliferation differences seen in the pooled Smart-seq2 analysis are confined to the primitive/HSC-labeled cluster (the biologically relevant compartment) or spread across more differentiated cell types in the plates.

## Conclusions / findings
Provides the per-cluster cycling/proliferation table used to hand-pick the LT-HSC-labeled cluster for the final matched-compartment comparison summarized in c_final.py.
