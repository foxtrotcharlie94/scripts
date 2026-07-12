# AgingCell_analysis.R

**Script:** `AgingCell_analysis.R`

## What it does
Downloads public dataset GSE207063 (Barman et al. 2022, Aging Cell) from GEO, filters to male Young vs Old classical-monocyte samples, demultiplexes via HTO, runs Seurat QC/clustering/UMAP, then pseudobulk-aggregates and runs edgeR quasi-likelihood DE (Old vs Young). Follows with GSEA (GO BP, KEGG, MSigDB Hallmarks) via clusterProfiler, producing bar/dot/enrichment-map/curve plots. The latter half cross-references this public aging signature against the researcher's own HB-vs-LB and Lenti+/- monocyte data, building custom "Old_up"/"Old_down" and "HB_up"/"HB_down" gene sets and running GSEA in both directions (aging signature scored on HB/LB data, and HB/LB signature scored on aging data), plus Lenti+ vs Lenti- within HB and within LB.

## Inputs
GSE207063 GEO supplementary files (10x matrices, HTO hashtags); own DE result objects your_mono (HB vs LB), lenti_hb/lenti_lb (Lenti+ vs Lenti- in HB/LB) presumably loaded earlier in an R session; res_df (aging DEGs).

## Outputs
GSE207063_seurat_male.rds, GSE207063_DEGs_male_OldVsYoung.csv, GSE207063_GSEA_GO_BP/KEGG/Hallmarks.csv, GSE207063_umap_male.png, GSE207063_volcano_male.png, various GSEA bar/dot/emap/curve PNGs, GSEA_OldvsYoung_on_HBvsLB.png, GSEA_OldvsYoung_on_LentiHB.png, GSEA_OldvsYoung_on_LentiLB.png, GSEA_HBvsLB_on_OldvsYoung.png.

## Conceptual aim / target
Test whether the transcriptional aging signature in monocytes overlaps with (or is recapitulated by) the clonal-hematopoiesis high-burden / lentiviral-marking transcriptional signatures — i.e., does CH-driven monocyte state resemble physiological aging.

## Note
Distinct analysis — only script using an external public GEO dataset; uniquely cross-links aging vs CH/lenti signatures via bidirectional GSEA.
