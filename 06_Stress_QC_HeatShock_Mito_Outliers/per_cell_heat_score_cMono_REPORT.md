# per_cell_heat_score_cMono.R

**Script:** `per_cell_heat_score_cMono.R`

## What it does
Identical pipeline to per_cell_heat_score_HSC.R but for Classical Monocytes: normalizes heat-shock UMIs per cell, z-scores per gene, averages into a module score, generates ridge/violin/CDF/rescue plots.

## Inputs
nUMIs_in_each_cells_alldata_res0.6_cluster0.csv and total_nUMIs_in_each_cell_alldata_res0.6_cluster0.csv from Classical_Monocytes (1).

## Outputs
heatshock_ridge/violin/cdf/HB4_rescue/rescue_allsamples.pdf/.png and heatshock_scores_per_cell.csv.

## Conceptual aim / target
Same as per_cell_heat_score_HSC.R but for monocytes.

## Note
Near-duplicate of per_cell_heat_score_HSC.R — same logic, different input population/paths.
