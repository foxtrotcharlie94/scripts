# BMDM library-composition & logFC tail-asymmetry diagnostics

**Script:** `bmdm_asymmetry_diagnostics.py`

## What it does
Checks whether the observed asymmetry in numbers of up- vs down-regulated genes across the four BMDM comparisons could be a normalization/library-composition artifact rather than real biology. Part 1 computes, per sample and per group, what fraction of the TMM-normalized CPM library is consumed by the top 1/5/10/20/50 genes. Part 2 plots the logFC distribution and a ranked-logFC "waterfall" per comparison to visualize tail one-sidedness.

## Inputs
edgeR TMM-normalized CPM CSVs and DEG bulk .txt files (columns incl. logFC, FDR, gene_name) for the four comparisons HB.neg_vs_LB.neg, HB.pos_vs_LB.pos, HB.pos_vs_HB.neg, LB.pos_vs_LB.neg, each in its own BMDM/<comparison>/ folder.

## Outputs
bmdm_composition_check.png/.csv (per-sample and per-group library dominance), bmdm_logFC_distribution.png (histograms + waterfall plots), plus console tables of composition and tail up/down splits.

## Conceptual aim / target
Rule out that the burden (HB vs LB) vs. lentiviral-construct (pos vs neg) transcriptional asymmetry is an artifact of a handful of very highly expressed transcripts dominating the normalized library, rather than a genuine broad shift in gene expression.

## Conclusions / findings
Diagnostic script — its output (composition curves, tail split counts) is used to confirm the asymmetry is not driven by single dominant genes; the specific per-comparison numbers are written to the CSV/console at runtime and were used as supporting evidence alongside the concordance and heatmap analyses below.
