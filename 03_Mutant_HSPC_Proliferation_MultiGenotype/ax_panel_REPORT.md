# Asxl1 10x panel-extraction step (per sample tag)

**Script:** `ax_panel.py`

**Status:** Pipeline step 1 of 2 (paired with ax_join.py)

## What it does
Per-sample preprocessing step for the Asxl1 GSE158184 10x dataset (4 HTO-co-captured samples: wk04_wt, wk04_mut, wk36_wt, wk36_mut). Loads one sample's matrix/features/barcodes, QC-filters, normalizes, scores cell cycle + proliferation-core, and either defines a shared highly-variable-gene panel ("define" mode, run on one sample) or extracts that same panel's expression for reuse (other samples), saving the panel subset matrix and per-cell metadata to a shared temp folder for later joining.

## Inputs
GSE158184 10x matrix.mtx.gz/features.tsv.gz/barcodes.tsv.gz per sample tag; sample tag, genotype, timepoint and mode ("define"/other) passed as command-line arguments.

## Outputs
/tmp/asxl1/panel.csv (shared gene panel, written once), /tmp/asxl1/<geno>_<timepoint>_X.npy and _obs.csv per sample.

## Conceptual aim / target
Extract a consistent, shared gene panel and per-cell cell-cycle/proliferation scores for each of the 4 Asxl1 samples individually, ahead of joining them for clustering and comparison in ax_join.py.

## Conclusions / findings
Preparation step only; per-sample cell counts and panel size are printed at runtime ("<geno> <tp> <mode> cells <n>"). Feeds directly into ax_join.py for the actual WT-vs-KO comparison.
