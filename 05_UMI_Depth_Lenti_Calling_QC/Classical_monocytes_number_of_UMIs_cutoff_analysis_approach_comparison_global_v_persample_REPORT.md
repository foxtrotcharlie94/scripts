# Classical_monocytes_number_of_UMIs_cutoff_analysis_approach_comparison_global_v_persample.R

**Script:** `Classical_monocytes_number_of_UMIs_cutoff_analysis_approach_comparison_global_v_persample.R`

## What it does
For LT-HSCs, fits a negative-binomial (fallback Poisson) regression of Lenti UMI counts on log(total nUMI) using two different training-set definitions: Approach A (global top 50% Lenti-frequency cells) and Approach B (top 50% per sample, then pooled), then derives the nUMI threshold at which P(detect >=1 Lenti UMI) = 80% for each approach and compares them.

## Inputs
LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt (gene x cell UMI matrix).

## Outputs
Lenti_threshold_comparison_LTHSCs.pdf/.png (two overlaid detection-probability curves with threshold annotations); printed threshold values.

## Conceptual aim / target
Decide/validate which cell-selection strategy (global vs per-sample ranking) should be used to define the lentiviral-detection UMI threshold for LT-HSCs, and quantify how much the resulting threshold differs between approaches.

## Note
Distinct analysis — only file focused on comparing threshold-derivation methodology itself rather than running downstream DE/GSEA.
