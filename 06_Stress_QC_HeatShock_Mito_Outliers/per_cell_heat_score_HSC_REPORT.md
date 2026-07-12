# per_cell_heat_score_HSC.R

**Script:** `per_cell_heat_score_HSC.R`

## What it does
For LT-HSCs, normalizes per-cell heat-shock gene UMI counts, z-scores each gene across cells, and averages into a "module score" per cell (hs_score); produces ridge/violin/CDF plots and an HB4-specific rescue analysis using z=0/z=1 thresholds.

## Inputs
nUMIs_in_each_cells_HSPCsubset_res0.3_cluster0.csv and total_nUMIs_in_each_cell_HSPCsubset_res0.3_cluster0.csv from LT-HSCs (1).

## Outputs
heatshock_ridge/violin/cdf/HB4_rescue/rescue_allsamples.pdf/.png and heatshock_scores_per_cell.csv.

## Conceptual aim / target
Quantify heat-shock/stress module activation per LT-HSC using a z-scored module score, comparing HB vs LB and flagging outlier samples.

## Note
Near-duplicate of per_cell_heat_score_cMono.R (same logic, different population); z-scored counterpart of heat_stress_score_per_cell_beta.R's raw %UMI fraction approach.
