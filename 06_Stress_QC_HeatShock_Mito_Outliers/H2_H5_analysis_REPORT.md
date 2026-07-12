# H2_H5_analysis.R

**Script:** `H2_H5_analysis.R`

## What it does
Investigates two candidate outlier samples (HB2, HB5) in LentiPos LT-HSC pseudobulk RNA-seq. Reads the LentiPosHB vs LentiPosLB pseudobulk CPM file, computes a PCA on log2(CPM+1) highlighting HB2/HB5, plots per-sample gene dropout, and computes % of total CPM contributed by proliferation genes vs mitochondrial OXPHOS genes to characterize HB5's transcript profile.

## Inputs
DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt (pseudobulk CPM matrix).

## Outputs
outlier_analysis_combined.pdf/.png (3-panel PCA + dropout + gene-set CPM-fraction figure), saved to "LT-HSCs (1)/LentiPosHB_vs_LentiPosLB".

## Conceptual aim / target
Determine whether HB2 and HB5 are technical/biological outliers in the LentiPos LT-HSC HB group, and whether HB5 shows an elevated proliferation/OXPHOS signature suggestive of technical aberration.

## Note
Nearly identical to the final "HB2_HB5_analysis" section embedded at the end of HB34_v_HB125_analysis.R — same code/logic/output filenames, essentially duplicated as a standalone file.
