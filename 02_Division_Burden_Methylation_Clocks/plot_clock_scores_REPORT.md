# Weighted clock scores by group, plus positive/negative-weight CpG split

**Script:** `plot_clock_scores.py`

## What it does
Plots the actual weighted Blood and WLMT clock scores (from apply_clocks.py's output) as grouped bar charts across the four biological groups (young/old baseline, young/old 10-cell reconstitution), and a second figure showing mean methylation split separately for age-up (positive-weight) vs. age-down (negative-weight) clock CpGs per group.

## Inputs
clock_cpg_analysis.csv (produced by apply_clocks.py).

## Outputs
clock_scores_by_group.png and clock_cpg_posneg_by_group.png; a printed summary table of clock scores by clock x group.

## Conceptual aim / target
Show, using the literature-calibrated clock scores (not just raw correlation), whether the 10-cell reconstitution group's predicted epigenetic age sits between the young and old baseline groups, consistent with accelerated aging from forced proliferation.

## Conclusions / findings
Provides the group-level bar-chart summary of clock score shifts referenced in the project's REPORT.md (see clock_correlation_analysis/ output folder) — old-recipient reconstitution group is explicitly labeled in the script as "not in Fig 2E", i.e. it is shown for context but excluded from the original paper's figure being validated.
