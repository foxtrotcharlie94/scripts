# Project overview — Human / Cross-Species GSEA

Extends the mouse clonal-hematopoiesis findings to human data: GSEA on human CD14+ monocyte
TET2-CHIP/CH differential expression tables (`DGE_tables/run_gsea.R` -> `plot_gsea_human.R`,
and the consolidated `all_genes_comparison.R`), and reciprocal signature-transfer testing whether
human TET2/aging/CH signatures are enriched in mouse HB/LB or LentiPos/Neg rankings (via
ortholog-mapping, `run_signature_gsea.R`, `pseudobulk_DE_GSEA/run_humanSignature_GSEA.R`, and
`LTHSC_pathway_PCA.R`'s six-dataset Jakobsen HSPC cross-species comparison).
