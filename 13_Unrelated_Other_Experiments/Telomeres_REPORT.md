# Telomeres.R

**Script:** `Telomeres.R`

## What it does
Fits a negative binomial regression of telomere length against mutation status in six CH driver genes (ASXL1, TET2, DNMT3A, TP53, JAK2, Multihit) plus a patient covariate, on a genotype-annotated human patient dataset; plots coefficient estimates with 95% CI as a forest/dot plot.

## Inputs
all_patients_genotype_collapsed_annotated_wide_v2.csv (human patient-level genotype + telomere length table).

## Outputs
Unsaved ggplot effect-size plot; no files written.

## Conceptual aim / target
Test whether specific clonal hematopoiesis driver mutations are associated with altered telomere length in human patients.

## Note
Unrelated to the mouse CH/lentiviral single-cell project — a human patient telomere-genotype association study.
