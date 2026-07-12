# Smart-seq2 build+cluster in one script (KMeans fallback)

**Script:** `c_ss2_full.py`

**Status:** Fallback rewrite of c1_build_ss2.py + c2_cluster_ss2.py

## What it does
A self-contained alternative to the c1+c2 two-step pipeline: builds, scores, scales/PCAs, and clusters the WT/Tet2/Dnmt3a Smart-seq2 data in a single run, using KMeans (6 clusters) instead of Leiden for the clustering step, then labels clusters by lineage marker and computes the same per-cluster-x-genotype cycling/proliferation table.

## Inputs
Same GSE124822 Smart-seq2 count files.

## Outputs
ss2_cluster_proliferation.csv, ss2_cluster_marker_scores.csv, ss2_obs_clustered.csv, fig4_ss2_pca.png.

## Conceptual aim / target
Provide a working alternative to the Leiden-based c1/c2 pipeline when the Leiden/igraph dependency was unavailable, using KMeans on PCA coordinates instead so the compartment-matched comparison could still be completed.

## Conclusions / findings
Reaches the same class of per-cluster cycling/proliferation result as c2_cluster_ss2.py via a different clustering algorithm — used as a fallback rather than a distinct conclusion.
