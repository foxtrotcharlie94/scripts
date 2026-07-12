# Jak2 (GSE227026) 10x panel-extraction step (per sample)

**Script:** `jak2_panel.py`

**Status:** Panel-extraction step for an upgraded dataset; no join script recovered

## What it does
Per-sample preprocessing step for a more replicated Jak2 dataset (GSE227026, part of the same multi-model preleukemic-mouse resource used for the upgraded Dnmt3a data). Reads a 10x .h5 file out of a per-sample tar archive, QC-filters, normalizes, scores cell cycle + proliferation-core, and defines/reuses a shared gene panel, labeling cells with a per-mouse ID.

## Inputs
<GSM>.tar containing one 10x .h5 file; GEO sample ID, genotype, replicate number and mode passed as command-line arguments.

## Outputs
/tmp/jak2/panel.csv (shared panel) and /tmp/jak2/<geno><rep>_X.npy/_obs.csv per sample.

## Conceptual aim / target
Prepare per-sample panels for the better-powered GSE227026 Jak2 comparison (3 WT vs 3 mutant replicates) identified as an upgrade over the original GSE263101 Jak2 data used in c_join.py.

## Conclusions / findings
Preparation step only; no corresponding join/comparison script for this specific GSE227026 Jak2 panel was recovered among the extracted files — the Jak2V617F conclusion reported in c_final.py instead comes from the GSE263101-based c_join.py pipeline.
