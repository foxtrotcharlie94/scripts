# Lenti_UMI_cutoff.R

**Script:** `Lenti_UMI_cutoff.R`

## What it does
Runs five independent statistical approaches to justify a UMI-depth cutoff for reliable Lenti+/Lenti- calling in LT-HSC and classical monocyte cells: binomial detection-probability model, empirical downsampling consistency check, Wilson confidence intervals per cell, a two-component Beta mixture model fit via EM, and an overlay summarizing all four approaches.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt and classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt.

## Outputs
detection_approach1_binomial.pdf through detection_approach5_summary.pdf, and lenti_detection_per_cell_scores.csv.

## Conceptual aim / target
Determine the minimum sequencing depth at which Lenti+ vs Lenti- classification becomes statistically reliable.

## Note
Companion/superset of Lenti_UMI_cutoff_50percent_pool.R (in name only — that file's actual content is unrelated, see its own entry).
