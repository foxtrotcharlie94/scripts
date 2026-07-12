# Smart-seq2 pipeline step 1: build merged, scored, PCA'd dataset

**Script:** `c1_build_ss2.py`

**Status:** Pipeline step 1 of 2 (paired with c2_cluster_ss2.py)

## What it does
First of a two-step compartment-matched Smart-seq2 pipeline. Merges the WT/Tet2/Dnmt3a plates into one AnnData object, normalizes, scores cell cycle + proliferation-core + six lineage-marker signatures (HSC/MPP/MkP/Ery/Mye/Lymph), computes highly-variable genes, scales and runs PCA, and saves the resulting object for the clustering step.

## Inputs
Same GSE124822 Smart-seq2 count files as s1_smartseq2.py.

## Outputs
ss2_built.h5ad (AnnData with scores + PCA embedding).

## Conceptual aim / target
Separate the (slower, one-time) data loading/scoring/PCA step from the (faster, iterated-on) clustering step, so clustering parameters can be tuned without re-processing the raw data each time.

## Conclusions / findings
Preparation step only — produces the input for c2_cluster_ss2.py, which generates the actual per-cluster comparison numbers.
