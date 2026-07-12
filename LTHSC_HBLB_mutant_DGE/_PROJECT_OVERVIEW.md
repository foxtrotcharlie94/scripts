# Project overview — HB/LB x Lenti+/- Differential Expression (Classical Monocytes & LT-HSCs)

The largest R sub-project: pseudobulk differential expression and volcano/summary plots testing
two crossed factors — clonal **burden** (HB = high-burden vs LB = low-burden clone size) and
**lentiviral-construct status** (LentiPos vs LentiNeg, i.e. whether a cell carries the clonal
marker) — in classical monocytes and LT-HSCs.

Structure: `plots_bmmono_*.R` / `plots_hsc_*.R` are single-comparison volcano plots (one script
per contrast, sharing an identical curated inflammation/aging or HSC-aging gene-label template);
`analysis1_HSC*.R` / `analysis1_CM.R` integrate all five comparisons per cell type into
overlap/effect-size/variance-partitioning summaries testing whether burden or genotype dominates;
`gene_and_pathway_analysis_beta_*.R` (7 files) run the more rigorous limma-voom/limma-CPM DGE +
5-database GSEA pipeline, each in two robustness arms (all samples vs excluding suspected
outlier animals HB3/HB4) and two cell types; `count_Genes_FDR.R`, `pathway_summary_table.R`,
`scatter_*.R` are LT-HSC-specific summary/scatter utilities (some duplicated across a
"LT-HSCs (1)" and "LT-HSCs (1) - Copy" folder pair — the Copy versions are byte-identical and
are included for completeness); `DEG_summary.R`, `concordance_analysis.R`,
`count_pathways_genes*.R` compare the "all samples" vs "excluding HB3/HB4" arms for robustness.
