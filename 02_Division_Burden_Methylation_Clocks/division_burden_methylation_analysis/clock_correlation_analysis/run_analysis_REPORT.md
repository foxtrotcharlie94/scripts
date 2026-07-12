# run_analysis.R

**Script:** `run_analysis.R`

## What it does
Reproduces (in R) a Python analysis behind Fig 2E of Gorelov et al. 2024, testing whether forced HSC proliferation mimics natural aging methylation. Extracts GSE44117 RRBS BED files, loads three published mouse DNAm clock CpG-coefficient tables (Blood/WLMT/YOMT), liftOvers clock coordinates mm10->mm9, computes per-CpG methylation deltas (old-vs-young; reconstituted-vs-young), and fits/plots correlations between the two deltas per clock.

## Inputs
GSE44117_RAW.tar (RRBS bed.gz files), elife-40675-supp3-v2.xlsx (clock CpG coefficients), mm10ToMm9.over.chain.gz.

## Outputs
cpg_deltas.csv, fit_stats.csv, overlay_all_clocks.png, and one <clock>_clock_correlation.png per clock.

## Conceptual aim / target
Test whether forced HSC division recapitulates natural aging-associated DNA methylation changes at clock CpGs.

## Note
The one script in the R codebase belonging to the division-burden/methylation-clock sub-project; reproduces a Python pipeline in R.
