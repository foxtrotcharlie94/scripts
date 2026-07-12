# H2B_analysis.R

**Script:** `H2B_analysis.R`

## What it does
Same pipeline as H2B_GSEA\H2B_analysis.R, but defines LR_up/LR_down using all genes with padj < 0.05 (split by direction) instead of a fixed top/bottom-200 cutoff, and annotates plots with the resulting gene-set sizes; also copies itself into OUTDIR for provenance.

## Inputs
DEgene2_analysis_Res_H2BGFP_Low_VS_High.csv and the same 5 LT-HSC pseudobulk DEG files under LT-HSCs (1)/<comparison>/.

## Outputs
Same output file set as H2B_GSEA\H2B_analysis.R, written to H2B_fdrcutoff/ instead.

## Conceptual aim / target
Same as H2B_GSEA\H2B_analysis.R, but using a statistically-driven (FDR-based) rather than fixed-size gene-set definition, as a robustness check.

## Note
Near-duplicate of H2B_GSEA\H2B_analysis.R — a matched pair testing sensitivity of the H2B GSEA results to gene-set-definition method.
