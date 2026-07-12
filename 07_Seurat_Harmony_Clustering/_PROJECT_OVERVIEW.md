# Project overview — Seurat/Harmony Clustering & UMAP Density Maps

Builds batch-corrected (Harmony) Seurat clustering/UMAP embeddings for LT-HSCs and classical
monocytes, then visualizes whether HB vs LB (and, in the "_Lenti" variants, Lenti+ vs Lenti-)
cells occupy different transcriptional neighborhoods via 2D kernel-density difference maps.
`cluser_composition_analysis.R` builds on the resulting Seurat objects to find cluster marker
genes and quantify per-cluster cell-frequency differences across the four HB/LB x Lenti+/-
groups. `LTHSC_pathway_PCA_Hallmark.R` adds a GSVA-pathway-level PCA with permutation/PERMANOVA
significance testing of group separation.
