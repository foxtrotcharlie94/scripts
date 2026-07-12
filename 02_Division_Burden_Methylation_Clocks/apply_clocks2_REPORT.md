# Apply YOMT + Liver clocks and compute per-CpG age/reconstitution deltas

**Script:** `apply_clocks2.py`

## What it does
Extends the clock analysis to two more published mouse clocks (YOMT and Liver), and — unlike apply_clocks.py — computes results per individual CpG rather than per sample: for each clock CpG, the mean methylation across young vs old baseline samples (age delta) and young baseline vs young 10-cell-reconstitution samples (reconstitution delta).

## Inputs
clocks2_mm9.pkl (YOMT + Liver clock CpG coordinates/weights) and the same 7 RRBS bed.gz files (2 young, 2 old, 3 reconstitution).

## Outputs
cpg_deltas_clocks2.csv — per-CpG clock, chr, pos, weight, delta_age_old_minus_young, delta_reconst_minus_young.

## Conceptual aim / target
Generate the per-CpG delta values needed to directly correlate natural aging change against forced-proliferation change at the same genomic sites, for two additional clocks beyond Blood/WLMT.

## Conclusions / findings
Sets up the data used later by cpg_delta_overlay3.py to add the YOMT clock into the combined 3-clock comparison (the Liver clock was computed but ultimately dropped from the downstream overlays).
