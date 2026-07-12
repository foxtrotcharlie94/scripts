# GSE44117 global methylation by group (Fig 2E setup)

**Script:** `analyze_gse44117.py`

## What it does
First-pass script on the reanalysis of GSE44117 (Beerman et al. 2013 mouse HSC RRBS). Assigns each sample to a biological group (fetal, young/old baseline, young/old 10-cell-reconstitution, 5-FU-treated) from its sample name, computes group-level mean +/- SD of global weighted CpG methylation, and plots a bar chart for the four groups relevant to Figure 2E of the paper being validated (young/old baseline vs. young/old 10-cell transplant reconstitution).

## Inputs
/tmp/gse44117/summary.csv — a pre-computed per-sample table (global weighted %methylation, n_cpgs, mean coverage) built earlier in the pipeline from the raw RRBS BED files.

## Outputs
summary_with_groups.csv, group_summary.csv, fig2e_groups_global_methylation.png.

## Conceptual aim / target
Establish a first, coarse (whole-genome average) comparison of methylation across natural aging (young vs old) and forced-proliferation reconstitution (10-cell transplant), as a stepping stone before applying actual epigenetic-clock CpG weights.

## Conclusions / findings
Explicitly flagged in the script's own plot title as "a proxy only — NOT the BS-WLMT/BS-Mouse-Blood clock scores in Fig 2E"; superseded by the clock-based scripts below for the actual hypothesis test.
