# Jak2 (GSE263101) 10x panel-extraction step used by c_join.py

**Script:** `t_panel.py`

**Status:** Pipeline step 1 of 2 (paired with c_join.py) — the version actually used for the final Jak2V617F result

## What it does
Per-genotype preprocessing step for the GSE263101 10x Jak2V617F dataset: extracts one genotype's matrix/features/barcodes from the dataset's tar archive, QC-filters, normalizes, scores cell cycle + proliferation-core, and defines ("define" mode) or reuses a shared HVG + marker gene panel, saving the panel-subset expression matrix and per-cell scores.

## Inputs
GSE263101_RAW.tar; working directory, genotype, file prefix and mode passed as command-line arguments.

## Outputs
<workdir>/tenx_panel.csv (shared panel, written once) and tenx_<genotype>_panel.npy/_obs.csv per genotype.

## Conceptual aim / target
Produce the exact per-genotype panel files consumed by c_join.py, i.e. this is the concrete panel-extraction step for the Jak2V617F result that made it into the project's final matched-compartment summary.

## Conclusions / findings
Preparation step only; per-genotype cell counts and panel size printed at runtime. Its output directly feeds the Jak2V617F row of c_final.py's MASTER_matched_compartment.csv via c_join.py.
