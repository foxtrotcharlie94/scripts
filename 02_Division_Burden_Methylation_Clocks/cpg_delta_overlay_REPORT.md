# Overlay figure: Blood + WLMT correlation on one panel (v1)

**Script:** `cpg_delta_overlay.py`

## What it does
Cosmetic/presentation iteration on cpg_delta_correlation.py: instead of two separate subplots (one per clock), plots both Blood and WLMT clock CpGs on a single shared axes so the two clocks' correlation lines can be compared directly, with per-clock stats annotated as text boxes.

## Inputs
cpg_deltas.csv (produced by cpg_delta_correlation.py).

## Outputs
delta_age_vs_delta_reconst_overlay.png.

## Conceptual aim / target
Make a single, compact figure comparing how strongly each clock's CpGs show the aging/reconstitution correlation, rather than requiring two side-by-side panels.

## Conclusions / findings
Purely a visualization iteration; the underlying r/slope/p values are identical to cpg_delta_correlation.py, just replotted. Superseded by cpg_delta_overlay3.py (adds a third clock) and overlay4.py (recolored).
