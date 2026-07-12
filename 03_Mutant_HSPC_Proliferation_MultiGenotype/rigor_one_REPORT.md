# Rigorous DESeq2 test for a single dataset (Asxl1 or Jak2)

**Script:** `rigor_one.py`

**Status:** Utility / debugging variant of rigor_bulk.py

## What it does
A single-dataset command-line variant of rigor_bulk.py: takes "asxl1" or "jak2" as an argument and runs the same DESeq2 + competitive Mann-Whitney / signed Wilcoxon proliferation gene-set test for just that one dataset, saving its own small CSV.

## Inputs
Same raw data files as rigor_bulk.py, but only the ones needed for the selected dataset; RES output directory and dataset name ("asxl1"/"jak2") passed as command-line arguments.

## Outputs
rigor_<dataset>.csv (single-dataset result row).

## Conceptual aim / target
Allow re-running or debugging the DESeq2 proliferation test for one dataset at a time without re-running the (slower) full 3-dataset rigor_bulk.py script.

## Conclusions / findings
Functionally a refactored subset of rigor_bulk.py; no additional conclusions beyond what rigor_bulk.py reports for the same two datasets.
