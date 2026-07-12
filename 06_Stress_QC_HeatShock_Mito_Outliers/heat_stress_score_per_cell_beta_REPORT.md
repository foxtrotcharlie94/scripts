# heat_stress_score_per_cell_beta.R

**Script:** `heat_stress_score_per_cell_beta.R`

## What it does
Loads per-cell heat-shock gene UMI counts and total UMI counts for Classical Monocytes then LT-HSCs, computes hs_fraction = HS-gene UMIs / total UMIs x 100 per cell, and produces ridge/violin/CDF plots plus an HB4-specific "rescue" analysis (fraction of HB4 cells below a p90 threshold derived from other samples).

## Inputs
nUMIs_in_each_cells_*_cluster0.csv and total_nUMIs_in_each_cell_*_cluster0.csv for Classical_Monocytes (1) and LT-HSCs (1).

## Outputs
heatshock_fraction_ridge/violin/cdf/HB4_rescue/rescue_allsamples.pdf/.png and heatshock_fraction_per_cell.csv.

## Conceptual aim / target
Determine whether high-burden clone samples (especially HB3/HB4) show elevated heat-shock/stress transcriptional signatures, and how many cells could be data-driven "rescued".

## Note
"Raw fraction" (%UMI) sibling of per_cell_heat_score_HSC.R/per_cell_heat_score_cMono.R, which compute a z-scored module score instead; nearly identical plot code/structure.
