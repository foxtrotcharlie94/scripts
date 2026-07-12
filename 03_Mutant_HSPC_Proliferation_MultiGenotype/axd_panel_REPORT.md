# Asxl1 10x panel extraction with in-matrix hashtag demultiplexing

**Script:** `axd_panel.py`

**Status:** Alternate ingestion path for ax_panel.py (adds HTO demultiplexing)

## What it does
A variant of ax_panel.py for an HTO-multiplexed version of the Asxl1 data: instead of relying on separate WT/KO sample files, reads sample hashtag (HTO) counts directly out of the combined feature matrix, assigns each cell to its biological replicate/genotype by argmax HTO signal (with a minimum-fraction/minimum-count doublet filter), then proceeds with the same QC, normalization, cell-cycle/proliferation scoring and panel-extraction logic as ax_panel.py.

## Inputs
10x matrix/features/barcodes for a given timepoint tag, where the features include per-replicate sample-hashtag rows matching pattern "wk<N>_(wt|mut)_<rep>"; timepoint tag and mode passed as command-line arguments.

## Outputs
/tmp/axd/panel.csv (shared panel) and /tmp/axd/<tag>_X.npy/_obs.csv (now including a "biorep" column identifying the specific mouse).

## Conceptual aim / target
Support an HTO-multiplexed ingestion path for Asxl1 data where genotype/replicate isn't already split into separate files, recovering per-mouse (not just per-genotype) resolution for more rigorous, replicate-aware statistics.

## Conclusions / findings
Preparation step; enables mouse-level (rather than only cell-level) pseudobulk statistics downstream, addressing pseudo-replication concerns in the plain ax_panel.py/ax_join.py cell-level comparisons.
