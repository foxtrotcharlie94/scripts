# Core hypothesis test: aging delta vs. reconstitution delta per CpG (Blood, WLMT)

**Script:** `cpg_delta_correlation.py`

## What it does
The central analysis of this project. For every Blood- and WLMT-clock CpG with full coverage across all 7 samples, computes the natural-aging methylation delta (old baseline minus young baseline) and the forced-proliferation delta (young 10-cell reconstitution minus young baseline), then correlates the two per clock (Pearson r, linear regression slope/intercept/p-value) in two separate scatter panels.

## Inputs
clocks_mm9.pkl and the 7 relevant RRBS bed.gz files (parses them directly rather than reusing apply_clocks.py's output).

## Outputs
cpg_deltas.csv (per-CpG deltas for both clocks) and delta_age_vs_delta_reconst_correlation.png (2-panel scatter with fit stats).

## Conceptual aim / target
Directly test the paper's central claim: if forced proliferation (10-cell transplant) drives epigenetic aging similarly to natural chronological aging, CpGs that gain/lose methylation with age should show a correlated gain/loss with reconstitution, at a slope informative about the relative "dose" of aging induced.

## Conclusions / findings
Reports, per clock, the Pearson r / regression slope / n CpGs / p-value printed at runtime (e.g. "Blood: slope=... r=... p=... n=..."); this is the key correlation figure referenced by the later overlay scripts, which is a positive correlation for the tested clocks based on how the downstream overlay scripts were iterated (higher slope/greater r for WLMT than Blood was consistently highlighted across script revisions).
