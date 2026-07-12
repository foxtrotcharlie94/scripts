# Jak2 (GSE263101) 10x: save full log-normalized matrix (all genes)

**Script:** `t_save.py`

**Status:** Utility export variant of t_panel.py (full matrix instead of panel subset)

## What it does
Near-duplicate of t_panel.py that, instead of saving only the shared gene-panel subset, saves the FULL log-normalized expression matrix (all genes, as a sparse .npz) plus the gene name list and per-cell cell-cycle/proliferation scores for one genotype.

## Inputs
GSE263101_RAW.tar; working directory, genotype and file prefix passed as command-line arguments.

## Outputs
<workdir>/tenx_<genotype>_lognorm.npz (full sparse matrix), tenx_<genotype>_obs.csv, tenx_varnames.csv.

## Conceptual aim / target
Preserve full-transcriptome access (not just the reduced marker/HVG panel) for a genotype's cells, in case a downstream step needed genes outside the shared panel (e.g. a follow-up differential-expression or signature check beyond what t_panel.py's subset supports).

## Conclusions / findings
Data-export utility; no comparison/statistics of its own. Complements t_panel.py's reduced-panel export for the Jak2 dataset.
