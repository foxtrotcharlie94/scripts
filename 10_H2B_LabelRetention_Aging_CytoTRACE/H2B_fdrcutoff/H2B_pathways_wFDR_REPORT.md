# H2B_pathways_wFDR.R

**Script:** `H2B_pathways_wFDR.R`

## What it does
A very large (4085-line), multi-stage pipeline built around the public GSE207063 dataset (Barman et al. 2022, Aging Cell). It downloads the GEO dataset, loads male 10x samples, runs Seurat QC/clustering/UMAP, pseudobulk-aggregates and runs edgeR DE (Old vs Young), then does GSEA (GO BP, KEGG, Hallmark) on the aging ranked list. It cross-references this aging signature against the researcher's own HB-vs-LB LT-HSC/monocyte Lenti+/Lenti- DEG tables via cross-dataset GSEA (forward and reverse), producing summary dotplots. The core block is duplicated/rerun ~5-6 times in sequence (apparent copy-paste iteration), and the file ends with a distinct final section testing a TET2-mutant signature (Jakobsen human HSPC/monocyte data) against the same HB-vs-LB contrasts.

## Inputs
GSE207063 GEO supplementary 10x files (downloaded via GEOquery); researcher's own pseudobulk DEG files for LT-HSC/monocyte HB vs LB and LentiPos vs LentiNeg (from "Classical_Monocytes (1)"); Jakobsen TET2 human DEG Excel files for the final section.

## Outputs
Writes to OUTDIR "AgingCell_analysis_wFDR" — DEG CSV, multiple GSEA CSVs/plots for GO BP/KEGG/Hallmarks, several "Summary_dotplot_crossdataset.png" files, and a separate OUTDIR_TET2 "TET2_signature_vs_HBvsLB" with TET2_GSEA_summary.csv and Summary_dotplot_TET2_on_HBvsLB.png.

## Conceptual aim / target
Determine whether physiological monocyte aging and/or a TET2-CH mutant signature recapitulate the lab's own HB-vs-LB transcriptional signatures in LT-HSCs and monocytes.

## Note
Internally repetitive — the core pipeline block is duplicated near-verbatim ~5-6 times within this single file (apparent iterative reruns); FDR-based sibling of H2B_GSEA/H2B_fdrcutoff H2B_analysis.R conceptually, but built around the GSE207063 aging dataset instead of the H2B label-retention dataset.
