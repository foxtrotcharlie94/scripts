# Asxl1 10x: self-contained single-script version of ax_panel.py + ax_join.py

**Script:** `asxl1_10x.py`

**Status:** Consolidated alternative to ax_panel.py + ax_join.py

## What it does
A self-contained, single-run alternative to the two-step ax_panel.py/ax_join.py pipeline: loads all 4 Asxl1 samples directly in one script, does cell-cycle scoring, cell-cycle-regressed PCA/KMeans clustering, lineage labeling, defines the LT-HSC+MPP1 compartment, and runs the same chi-square/Mann-Whitney/paired-pseudobulk comparisons as ax_join.py.

## Inputs
GSE158184 10x matrix/features/barcodes files for all 4 sample tags (wk04_wt, wk04_mut, wk36_wt, wk36_mut), located anywhere under the project's data folder (found via recursive glob).

## Outputs
asxl1_10x_obs.csv; console cluster composition, per-timepoint/pooled comparison stats, and paired-pseudobulk test result.

## Conceptual aim / target
Provide a single, easier-to-rerun script for the Asxl1 compartment-matched comparison, avoiding the need to coordinate 4 separate ax_panel.py invocations plus a join step.

## Conclusions / findings
Reaches the same class of WT-vs-KO LT-HSC+MPP1 result as ax_join.py in one script; effectively a consolidated rewrite of the ax_panel/ax_join pair for convenience.
