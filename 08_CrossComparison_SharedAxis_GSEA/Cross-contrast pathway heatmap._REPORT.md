# Cross-contrast pathway heatmap..R

**Script:** `Cross-contrast pathway heatmap..R`

## What it does
Reads pre-computed GSEA CSVs for 4 contrasts (LentiNeg_HBvsLB, LentiPos_HBvsLB, HB_PosVsNeg, LB_PosVsNeg) across 2 gene-set collections (Hallmark, GO_BP) and 2 rankings (logFC, logFCxP), builds the union of pathways significant (FDR<0.05) in at least one contrast, ranks them by mean(NES x -log10FDR), splits into up/down blocks, and draws a rows=pathways x cols=contrasts NES heatmap with significance stars.

## Inputs
GSEA_<collection>_<comparison>_rankBy_<ranking>.csv files under pop_dir = pseudobulk_DE_GSEA/ClassicalMonocytes/<comparison>/.

## Outputs
heatmap_table_union_<coll>_rankBy_<rk>.csv and heatmap_union_signature_<coll>_rankBy_<rk>.png per collection x ranking combination, written to pop_dir/_shared_axis/.

## Conceptual aim / target
Visualize which Hallmark/GO_BP pathways are shared or distinct across the four burden (HB/LB) vs lentiviral-marking (Pos/Neg) contrasts in classical monocytes.

## Note
Distinct analysis — a lightweight downstream visualization consumer of GSEA CSVs (contrasts with "plot_gsea_human.R" which operates on the human CD14 monocyte cohort instead).
