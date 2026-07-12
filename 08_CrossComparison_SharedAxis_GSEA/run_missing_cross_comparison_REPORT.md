# run_missing_cross_comparison.R

**Script:** `run_missing_cross_comparison.R`

## What it does
Driver/orchestration script that scans all pseudobulk_DE_GSEA/<population> folders, identifies which have the 4 required GSEA comparison subfolders but are missing a _cross_comparison output, and for each missing one sources an external DE_GSEA_analysis.R script to generate it.

## Inputs
Directory structure under pseudobulk_DE_GSEA/ and the external DE_GSEA_analysis.R script.

## Outputs
Triggers creation of _cross_comparison/ subfolders (same kind of output as gsea_cross_comparison.R) for each population missing them.

## Conceptual aim / target
Batch-automate re-running the cross-comparison GSEA analysis across every cell-population folder not yet processed.

## Note
Orchestration wrapper around gsea_cross_comparison.R's logic rather than a duplicate analysis itself.
