# Generalized, stricter compartment-matched join/comparison script

**Script:** `strict_join.py`

**Status:** General-purpose, stricter successor to ax_join.py / c_join.py

## What it does
A more general and stricter version of the *_join.py scripts: works for any genotype/dataset whose panel files follow the shared /<TMP>/<lib>_X.npy + _obs.csv convention, and for either "biorep" or "mouse" as the replicate-grouping unit. Uses finer KMeans clustering (20 clusters vs. 12 elsewhere), a stricter LT-HSC+MPP1 definition (HSC z-score high AND every lineage-signature z-score low, rather than just "HSC is the top-scoring signature"), drops low-cell-count units (<30 cells), and runs (optionally per-timepoint, then pooled) t-tests on cycling% and proliferation score between WT and any of KO/Dnmt3a/Jak2V617F.

## Inputs
Panel .npy/_obs.csv files for a given library set (path passed as command-line TMP argument) plus a matching panel.csv, in the same format produced by ax_panel.py / t2_panel.py / d_panel.py / jak2_panel.py / r878_panel.py.

## Outputs
<TMP>/strict_obs.csv (per-cell cluster/compartment labels) and <TMP>/strict_pb.csv (per-unit pseudobulk cycling/proliferation); console cluster table and WT-vs-mutant test results.

## Conceptual aim / target
Provide one reusable, more rigorous compartment-definition and statistical-testing script that can be pointed at any of the project's per-genotype panel outputs, rather than maintaining a separate bespoke *_join.py per mutant.

## Conclusions / findings
Intended as the project's more defensible, general-purpose replacement for the earlier per-genotype ax_join.py/c_join.py scripts (stricter compartment definition, replicate-level rather than cell-level testing, minimum-n filtering) — it is unclear from the recovered files whether it was ultimately re-run against every genotype's panel output or only used for a subset.
