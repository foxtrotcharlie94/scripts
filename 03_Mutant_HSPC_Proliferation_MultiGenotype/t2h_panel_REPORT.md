# Tet2 (GSE209994) 10x panel-extraction step — near-duplicate of t2_panel.py

**Script:** `t2h_panel.py`

**Status:** Near-duplicate of t2_panel.py — likely a debugging checkpoint copy

## What it does
Line-for-line functionally identical to t2_panel.py (same HTO demultiplexing, QC, scoring and panel-extraction logic), writing to a differently-named temp directory (/tmp/tet2h instead of /tmp/tet2) — consistent with a checkpoint/backup copy made during iterative debugging rather than a deliberate methodological variant.

## Inputs
Same GSE209994 lane files as t2_panel.py.

## Outputs
/tmp/tet2h/panel.csv and /tmp/tet2h/<group>_<lane>_X.npy/_obs.csv.

## Conceptual aim / target
Unclear distinct aim from t2_panel.py — most likely kept as a working checkpoint copy while iterating on the HTO-demultiplexing logic.

## Conclusions / findings
No distinct conclusions beyond what t2_panel.py produces; recommend treating t2_panel.py as the canonical version and t2h_panel.py as a redundant duplicate unless a specific difference in downstream files is identified.
