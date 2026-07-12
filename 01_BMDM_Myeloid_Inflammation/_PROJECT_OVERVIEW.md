# Project overview — BMDM Myeloid-Inflammation / Clonal-Hematopoiesis Burden

**Biological question:** in bone-marrow-derived macrophages (BMDM) carrying a
lentiviral clonal-hematopoiesis (CH) construct, does clone "burden" (HB =
high-burden vs LB = low-burden) and/or the lentiviral construct itself
(lentipos vs lentineg, donor-matched) drive a myeloid-inflammation
transcriptional program?

**Design:** 4 pairwise bulk RNA-seq comparisons (HB.neg_vs_LB.neg,
HB.pos_vs_LB.pos, HB.pos_vs_HB.neg, LB.pos_vs_LB.neg), each already processed
through edgeR (TMM-normalized CPM + DEG tables) in R. The 3 Python scripts in
this folder are downstream diagnostics and visualizations built on top of
those edgeR outputs.

**Scripts (read in this order):**
1. `bmdm_asymmetry_diagnostics.py` — rules out a normalization artifact.
2. `bmdm_axis_concordance.py` — tests whether burden and lenti effects share a
   transcriptional axis.
3. `ch_myeloid_heatmap.py` — pinpoints which specific inflammatory genes are
   significantly affected by each axis.

**Note:** the actual differential-expression modeling (edgeR/TMM) was done in
R (`BMDM_analysis_all.R`, `BMDM_volcano_DEcounts.R`, not included here since
this portfolio covers Python scripts only); these three Python scripts operate
on that R output.
