# LTHSC_pathway_PCA_Hallmark.R

**Script:** `LTHSC_pathway_PCA_Hallmark.R`

## What it does
Three-part script. Part 1 loads the classical-monocyte gene x cell count matrix (despite the folder name), classifies Lenti+/Lenti- per cell, filters genes detected in >=5% of cells per sample, pseudobulk-aggregates per sample x Lenti-status, TMM-normalizes to logCPM, computes GSVA pathway scores against MSigDB Hallmark, runs PCA on pathway scores, and plots samples plus top PC1/PC2 pathway loadings. Part 2 runs a permutation test (10,000 shuffles) on HB-vs-LB centroid distance in PCA space. Part 3 runs a paired PERMANOVA (vegan::adonis2) on PC1-4.

## Inputs
classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt (Part 1); pca_coords.csv from a prior PCA run (Parts 2-3).

## Outputs
PCA_samples.png, Loadings_PC1/2.png, gsva_scores.rds, pca_res.rds, gsva_scores_samplesXpathways.csv, pca_coords.csv (Part 1); PCA_separation.png, Permutation_distribution.png, stats_results.txt (Part 2); permanova_paired_results.txt (Part 3).

## Conceptual aim / target
Test whether pathway-level (Hallmark GSVA) transcriptional states separate HB from LB clone-burden groups in a statistically rigorous way (permutation test, PERMANOVA), beyond visual PCA clustering.

## Note
Distinct analysis combining GSVA/PCA with formal permutation/PERMANOVA testing; internal inconsistency — despite referencing LT-HSC in its name, Part 1 actually reads the classical_monocytes count matrix.
