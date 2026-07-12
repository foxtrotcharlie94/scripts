# Tet2 (GSE209994) 10x panel-extraction step with HTO demultiplexing

**Script:** `t2_panel.py`

**Status:** Pipeline step feeding the tet2_analysis/ Harmony clustering results (see README.md in that folder)

## What it does
Per-sample-lane preprocessing step for the Tet2 GSE209994 10x dataset (HTO-multiplexed, 2 biological replicates per lane, Untreated + IL-1b arms). Loads a lane's matrix/features/barcodes, demultiplexes hashtag (HTO) reads to assign each cell to one of its 2 replicate mice, filters low-confidence HTO assignments and low-quality cells, normalizes, scores cell cycle + proliferation-core, and defines/reuses a shared gene panel.

## Inputs
GSE209994_<group>_<lane>_filtered_feature_bc_matrix files (matrix/features/barcodes, features including 2 Antibody-Capture HTO rows); group, lane, genotype, condition and mode passed as command-line arguments.

## Outputs
/tmp/tet2/panel.csv (shared panel) and /tmp/tet2/<group>_<lane>_X.npy/_obs.csv, with a per-cell "mouse" ID from HTO demultiplexing.

## Conceptual aim / target
Recover per-mouse (not just per-genotype) resolution for the Tet2 baseline (Untreated) and IL-1b-stress comparisons, matching the "powered, batch-free" dataset identified in BETTER_DATASETS_Tet2_Dnmt3a.md as an upgrade over the original Smart-seq2 Tet2 data.

## Conclusions / findings
Preparation step; per-lane cell counts and per-mouse composition printed at runtime. This is the version whose results are described in the companion tet2_analysis/README.md as showing a small but consistent cross-cluster proliferation increase in Tet2-KO (sign-test p=0.007).
