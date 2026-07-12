# HB34_v_HB125_analysis.R

**Script:** `HB34_v_HB125_analysis.R`

## What it does
A large multi-section script comparing outlier samples HB3/HB4 against HB1/HB2/HB5 across four pseudobulk datasets (LT-HSC and Classical Monocyte, LentiPos/LentiNeg). Sections: PCA; limma DEG + fgsea GSEA volcano plots (Hallmark); individual-gene volcano plots; ORA on significant DEGs; heat/stress-vs-RNA-handling-artifact colored volcano plot; pheatmap of upregulated heat/stress genes; a "heat & stress score" bar chart; and an appended copy of the HB2/HB5 outlier analysis.

## Inputs
Four pseudobulk CPM DEG text files (LT-HSCs and Classical Monocytes, LentiNeg and LentiPos HB-vs-LB).

## Outputs
pca_pseudobulk.pdf/png, gsea_volcano_HB34_vs_HB125.pdf/png, gene_volcano_HB34_vs_HB125.pdf/png, ora_pathway_HB34_vs_HB125.pdf/png, heatmaps_heat_stress_up.pdf/png, heat_stress_score_all.pdf/png, heat_stress_score_hsp.pdf/png, plus outlier_analysis_combined.pdf/png.

## Conceptual aim / target
Characterize whether HB3 and HB4 are technical/biological outliers relative to HB1/2/5, and whether their divergent transcriptome reflects a heat-shock/stress-response or RNA-handling artifact rather than true biology.

## Note
Effectively 7 stacked mini-analyses in one file; its final section duplicates H2_H5_analysis.R almost verbatim.
