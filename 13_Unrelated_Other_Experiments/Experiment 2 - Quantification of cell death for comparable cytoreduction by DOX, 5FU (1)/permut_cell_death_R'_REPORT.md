# permut_cell_death_R'.R

**Script:** `permut_cell_death_R'.R`

## What it does
Runs a simple two-sample permutation test (100,000 shuffles) comparing sums of a per-mouse measurement between a DOX-treated group and a 5FU-treated group, separately for two cell populations (LT-HSC, and a "Lin-negative" population), computing an empirical two-sided p-value and plotting a histogram of the permuted null distribution with the observed difference marked.

## Inputs
Hardcoded numeric vectors (group_A = DOX, group_B = 5FU) pasted directly into the script; no external file read.

## Outputs
Base-R histogram plots (not saved to file); printed observed difference and p-value to console.

## Conceptual aim / target
Test whether DOX and 5-FU chemotherapy produce comparably-sized cytoreduction in LT-HSC and Lin-negative compartments, to justify treating them as comparable myeloablative conditions.

## Note
Distinct analysis — unrelated to the GSEA/DGE monocyte-HSC scripts; a small self-contained wet-lab statistics script from a separate cytoreduction-validation experiment.
