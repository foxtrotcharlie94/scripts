# BMDM_analysis_all.R

**Script:** `BMDM_analysis_all.R`

## What it does
A comprehensive diagnostic script over four pre-computed BMDM edgeR comparisons (HB.neg_vs_LB.neg, HB.pos_vs_LB.pos, HB.pos_vs_HB.neg, LB.pos_vs_LB.neg), reading each comparison's DEG table and TMM-normalized CPM matrix. It runs six analysis sections: (1) DGE up/down asymmetry at FDR thresholds, (2) logFC distribution/tail-asymmetry histograms and waterfall plots, (3) library-composition/top-N-gene dominance checks, (4) logFC vs logCPM intensity-dependence, (5) cross-axis concordance scatterplots comparing "burden" (HB vs LB) effect vs "lenti" (pos vs neg) effect for all genes and for high-movers, and (6) a curated CH/inflammation gene panel heatmap with paired/unpaired BH-corrected t-tests across the 4 groups.

## Inputs
BMDM/<comp>/DEG_<comp>_bulk.txt and BMDM/<comp>/edgeR_<comp>_TMM_normalized_CPM_wOffset.csv for each of the 4 comparisons.

## Outputs
Files under BMDM/axis_analysis/: 1_DGE_asymmetry_summary.csv/.png, 2_logFC_tail_split.csv, 2_logFC_distribution.png, 3_composition_per_sample.csv, 3_composition_check.png, 3b_abundant_gene_leverage.csv, 4_intensity_dependence.csv/.png, 5_axis_logFC.csv, 5a_concordance_allgenes.png, 5b_concordance_movers.png, plus a CH-gene heatmap (expression + effect-size panels).

## Conceptual aim / target
Rigorously characterize BMDM transcriptional response to clonal burden (HB/LB) versus lentiviral-marking status (pos/neg), check for technical artifacts (library composition, intensity bias), and assess whether the two axes (burden vs lenti) drive concordant or independent inflammatory gene programs.

## Note
Distinct, most extensive/QC-oriented BMDM script; complements BMDM_volcano_DEcounts.R (which only does volcano+bar plots) by adding deep diagnostic/concordance analysis on the same 4 comparisons.
