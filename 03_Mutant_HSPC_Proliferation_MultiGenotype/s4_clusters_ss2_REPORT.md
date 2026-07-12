# Smart-seq2: add Leiden clustering + lineage annotation to WT/Tet2/Dnmt3a

**Script:** `s4_clusters_ss2.py`

**Status:** Superseded by c1_build_ss2.py + c2_cluster_ss2.py

## What it does
Extends the pooled Smart-seq2 analysis by adding Leiden clustering, UMAP embedding, and lineage-marker-based cluster annotation (HSC/MPP/MkP/Ery/Mye/Lymph signature scores), so that %cycling and proliferation score can be reported per cell-type cluster x genotype instead of only pooled across all cells.

## Inputs
Same GSE124822 Smart-seq2 count files as s1_smartseq2.py.

## Outputs
ss2_cluster_proliferation.csv, ss2_cluster_marker_scores.csv, fig4_ss2_umap.png.

## Conceptual aim / target
Test whether any pooled-cell proliferation differences (from s1_smartseq2.py) are actually located in the primitive HSC-like cluster or are driven by more differentiated contaminating cells in the plates.

## Conclusions / findings
An intermediate clustering step; its Leiden-based approach was replaced in c1_build_ss2.py/c2_cluster_ss2.py by a more complete two-step build+cluster pipeline, and in c_ss2_full.py by a KMeans-based self-contained rewrite (avoiding a Leiden/igraph dependency issue).
