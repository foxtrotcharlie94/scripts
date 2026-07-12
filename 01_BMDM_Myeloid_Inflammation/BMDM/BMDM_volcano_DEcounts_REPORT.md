# BMDM_volcano_DEcounts.R

**Script:** `BMDM_volcano_DEcounts.R`

## What it does
Reads the same four BMDM DEG tables as BMDM_analysis_all.R and generates volcano plots for each comparison under three significance criteria (unadjusted P<0.05, FDR<0.25, FDR<0.05), saving both individual and combined faceted figures, with top genes labeled via ggrepel. Also produces a bar chart / CSV of total up/down DE gene counts per comparison at unadjusted P<0.05.

## Inputs
BMDM/<comp>/DEG_<comp>_bulk.txt for the 4 comparisons (HB.neg_vs_LB.neg, HB.pos_vs_LB.pos, HB.pos_vs_HB.neg, LB.pos_vs_LB.neg).

## Outputs
BMDM/volcano_DEcounts/volcano_<comp>_<crit>.png (individual), volcano_all_<crit>.png (combined faceted, x3 criteria), DE_counts_P0.05.csv, and an up/down bar plot.

## Conceptual aim / target
Visualize and quantify differential expression magnitude/direction in BMDMs across the burden and lenti-status axes at multiple significance stringencies.

## Note
Companion/simpler variant of BMDM_analysis_all.R — shares data-reading logic and the 4 comparisons, but restricted to volcano plots + DE counts (no composition/concordance/heatmap sections).
