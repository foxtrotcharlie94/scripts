# Project overview — Cross-Comparison / Shared-Axis GSEA

Addresses a central methodological question of the whole HB/LB x Lenti+/- study: how much does
the "burden" (HB vs LB) transcriptional signature overlap with the "lentiviral marking" (Pos vs
Neg) signature? `pseudobulk_DE_GSEA.R` is the core DE/GSEA engine producing the per-comparison
CSVs that `gsea_cross_comparison.R`, `shared_axis_heatmap.R`, `shared_axis_pathway.R`, and
`DE_GSEA_analysis.R` all consume, using different statistical framings (NES correlation and
meta-GSEA enrichment, static union heatmaps, and total-least-squares "dose" regression,
respectively) to test and quantify shared vs independent signal between the two axes.
`run_missing_cross_comparison.R` batch-orchestrates this analysis across cell populations.
