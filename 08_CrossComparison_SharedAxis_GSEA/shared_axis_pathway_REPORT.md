# shared_axis_pathway.R

**Script:** `shared_axis_pathway.R`

## What it does
Merges GSEA results for the two "burden" and two "lenti" contrasts, filters to pathways significant in at least one and consistent (FDR<0.10) across all four, fits a total-least-squares-through-origin regression of lenti NES on burden NES (with bootstrap CI) to estimate a "dose", and decomposes each contrast's pathway-response vector into the fraction aligned with the shared burden axis.

## Inputs
Same GSEA_<collection>_<comparison>_rankBy_<ranking>.csv files as shared_axis_heatmap.R.

## Outputs
shared_axis_<coll>_rankBy_<rk>.png (scatter+decomposition), shared_axis_byCondition_<coll>_rankBy_<rk>.png, shared_axis_dose.csv, shared_axis_decomposition.csv.

## Conceptual aim / target
Test/quantify the hypothesis that the LentiPos-vs-LentiNeg pathway signature is a smaller, same-direction move along the same axis as the HB-vs-LB burden signature.

## Note
Companion to shared_axis_heatmap.R; conceptually the same "shared axis" hypothesis as gsea_cross_comparison.R's captured-fraction analysis, implemented independently via TLS-through-origin regression.
