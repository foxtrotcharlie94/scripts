# Project overview — Barcode Lineage-Tracing & Progenitor-to-Tissue Distance

A distinct sub-project (different lentiviral barcode dataset, DARLIN-style lineage tracing, not
the HB/LB burden study) asking which bone-marrow progenitor population (LSK, MPP4, CLP, GMP,
across several bones) is the closest clonal source for various downstream tissue destinations
(thymus, blood T cells, blood/heart granulocytes). Utility scripts (`files/combine_barcodes.R`,
`combine_and_analyze_CT.R`, `barcode_visualizations.R`) build and visualize the barcode-sharing
tables; `distance_analyses_all_metrics.R` and `pooled_perm_tests.R` generalize the
distance/permutation-test logic to multiple metrics and datasets; `combined_regions_analysis*.R`
and `separate_regions_analysis.R` are earlier, less-generalized versions of the same pipeline;
`pooled_bootstrap_CLP_vs_GMP_destfixed.R` / `pooled_perm_CLP_vs_GMP_multi_distance.R` focus
specifically on the CLP-vs-GMP-to-T-lineage question with two complementary inference methods.
`lenti_overlaps_with_mouse_genome.R` and `neutral_region_pipeline.R` are related genome-annotation
utilities (vector-genome homology check; neutral safe-harbor region selection for a separate
liver mutagenesis assay) grouped here as adjacent barcode/vector-design tooling.
