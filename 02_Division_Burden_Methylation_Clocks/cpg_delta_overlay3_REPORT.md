# Overlay figure: Blood + WLMT + YOMT correlation (v2, 3 clocks)

**Script:** `cpg_delta_overlay3.py`

## What it does
Extends cpg_delta_overlay.py by merging in the YOMT clock's per-CpG deltas (from apply_clocks2.py's output, dropping the Liver clock) so all three clocks appear on one overlaid aging-delta vs. reconstitution-delta scatter with per-clock regression lines and stats.

## Inputs
cpg_deltas.csv (Blood, WLMT) and cpg_deltas_clocks2.csv (YOMT, Liver — Liver rows filtered out).

## Outputs
cpg_deltas_all3.csv (merged 3-clock table) and delta_age_vs_delta_reconst_overlay3.png.

## Conceptual aim / target
Broaden the epigenetic-clock validation from 2 to 3 independently published mouse clocks, checking whether the aging/reconstitution correlation is a general clock property rather than specific to one clock's CpG set.

## Conclusions / findings
Visualization/data-merging iteration; adds YOMT alongside Blood/WLMT with its own r/slope/p, again not hardcoded but printed and drawn at runtime. Superseded cosmetically by cpg_delta_overlay4.py.
