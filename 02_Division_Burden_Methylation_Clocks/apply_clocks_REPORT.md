# Apply Blood + WLMT epigenetic clocks to GSE44117 RRBS data

**Script:** `apply_clocks.py`

## What it does
Scans each sample's gzipped RRBS BED file for the genomic coordinates of the Blood clock (Petkovich et al., 90 CpGs) and WLMT clock (Meer et al., 435 CpGs), matching on lifted (mm10->mm9) chromosome/position with a +/-1bp offset tolerance, and computes per-sample weighted methylation clock scores plus coverage statistics.

## Inputs
clocks_mm9.pkl (pickled list of (chr, mm9 position, weight) tuples per clock, produced by an earlier liftOver step) and 9 per-sample RRBS bed.gz files (young/old baseline + 10-cell reconstitution groups).

## Outputs
clock_cpg_analysis.csv — one row per sample per clock, with coverage %, mean methylation split by positive/negative CpG weight, and the normalized weighted clock score.

## Conceptual aim / target
Compute the actual, literature-defined epigenetic-age clock scores (not just global methylation) for each sample, as the basis for testing whether forced proliferation shifts predicted epigenetic age.

## Conclusions / findings
Produces the per-sample clock scores later visualized in plot_clock_scores.py; coverage was good enough for both clocks to compute scores for all 9 samples (see script's own printed "done <sample>" log).
