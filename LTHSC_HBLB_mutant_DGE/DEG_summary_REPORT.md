# DEG_summary.R

**Script:** `DEG_summary.R`

## What it does
Compares limma pseudobulk DE results between "all samples" and "excluding HB3/HB4" analysis runs for LT-HSCs and Classical Monocytes, across 4 comparisons. For each comparison it produces a DEG-count bar chart (FDR<0.2), a logFC concordance scatter, and an FDR concordance scatter (labeling top-discordant genes), then assembles these into combined multi-panel PDFs per dataset.

## Inputs
limma log-CPM result CSVs (*_limma_*results.csv) under allSamples_analysis/ and excluding_H3_H4/ subfolders for both LT-HSCs and Classical Monocytes.

## Outputs
DEG_comparison_counts.pdf/.png, DEG_comparison_logFC.pdf/.png, DEG_comparison_FDR.pdf/.png, written to each dataset's summary_plots/ folder.

## Conceptual aim / target
Assess how sensitive the DE results are to excluding two specific samples (HB3, HB4) that may be outliers or batch-affected.

## Note
Distinct analysis — a sample-exclusion sensitivity check.
