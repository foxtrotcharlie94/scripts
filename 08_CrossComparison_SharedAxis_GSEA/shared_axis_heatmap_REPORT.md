# shared_axis_heatmap.R

**Script:** `shared_axis_heatmap.R`

## What it does
Reads saved GSEA CSVs for four HB/LB/Lenti contrasts within a given cell population, takes the union of pathways significant (FDR<0.05) in at least one, orders them by cross-contrast aggregate score, and draws a pathway x contrast heatmap with significance stars and a visual separator between "burden" and "lenti" columns.

## Inputs
GSEA_<collection>_<comparison>_rankBy_<ranking>.csv files under pop_dir/<comparison>/, for Hallmark and GO_BP, both rankings.

## Outputs
heatmap_table_union_<coll>_rankBy_<rk>.csv and heatmap_union_signature_<coll>_rankBy_<rk>.png.

## Conceptual aim / target
Visualize, across all four HB/LB x LentiPos/Neg contrasts simultaneously, which pathways are significantly changed and in what direction.

## Note
Sibling of shared_axis_pathway.R (same pop_dir/contrasts/output folder) — this is the heatmap view, the other is the scatter/decomposition ("dose") view.
