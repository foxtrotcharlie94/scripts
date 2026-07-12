# Burden-axis vs. lenti-axis concordance analysis

**Script:** `bmdm_axis_concordance.py`

## What it does
Tests whether two independent effect axes move together at the gene level: the "burden" axis (HB vs LB, averaged across the neg- and pos-background comparisons) and the "lenti" axis (lentipos vs lentineg, averaged across the HB- and LB-background comparisons). Produces a hexbin density scatter with Pearson r, Spearman rho, R-squared and a regression line for all genes, then a second scatter restricted to genes moving appreciably on both axes (|log2FC|>0.5), with per-quadrant percentages.

## Inputs
DEG bulk .txt files (logFC, gene_name) for the same four comparisons used by bmdm_asymmetry_diagnostics.py.

## Outputs
bmdm_concordance_allgenes.png, bmdm_concordance_movers.png; printed fit statistics (r, rho, R^2, slope) and quadrant percentages.

## Conceptual aim / target
Determine whether the CH "burden" (mutant clone size) and the lentiviral marking construct itself act through a shared transcriptional program (same genes moving together on both axes) or represent largely independent/orthogonal effects on BMDM gene expression.

## Conclusions / findings
Quantifies the fraction of "mover" genes that change in the same direction on both axes (printed "concordant = X%" plus per-quadrant counts) — the actual percentage is data-dependent and reported at runtime/in the saved PNGs rather than hardcoded.
