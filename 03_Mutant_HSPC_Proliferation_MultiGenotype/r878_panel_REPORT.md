# Dnmt3a-R878H 10x panel-extraction step with full-lineage marker panel

**Script:** `r878_panel.py`

**Status:** Panel-extraction step; no join/comparison script recovered

## What it does
Per-sample preprocessing step for a Dnmt3a-R878H point-mutant 10x dataset, using a substantially larger marker gene list than the other panel scripts — extending beyond HSC/MPP/Ery/Mye/Lymph into mature neutrophil, monocyte/macrophage, dendritic-cell, mast/basophil, B-cell and T/NK-cell markers — suggesting this dataset captures a broader, more mature cell-type range than the LSK-focused samples elsewhere in the project. QC-filters, normalizes, scores cell cycle + proliferation-core, and defines/reuses a shared panel.

## Inputs
<GSM>.tar (10x matrix/features/barcodes); GEO sample ID, genotype, mouse ID and mode passed as command-line arguments.

## Outputs
/tmp/r878/panel.csv (shared panel) and /tmp/r878/<mouse>_X.npy/_obs.csv per sample.

## Conceptual aim / target
Prepare per-sample panels for a Dnmt3a-R878H point-mutation comparison spanning the full myeloid/lymphoid maturation hierarchy, not just primitive HSPCs.

## Conclusions / findings
Preparation step only; no corresponding join/comparison script for R878H was recovered among the extracted files.
