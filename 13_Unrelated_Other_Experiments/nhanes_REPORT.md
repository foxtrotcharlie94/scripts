# nhanes.R

**Script:** `nhanes.R`

## What it does
Downloads multiple NHANES cycles (2005-2023) of CBC and demographic data via nhanesA, merges them, applies extensive exclusion criteria (cancer, anemia, kidney/liver disease, HIV, pregnancy, elevated CRP, smoking, obesity, alcohol, diabetes), then plots each CBC parameter against age, stratified by sex, using survey-weighted LOESS and quantile-regression IQR ribbons.

## Inputs
Downloaded directly from NHANES via nhanesA (CBC, DEMO, MCQ, KIQ, RHQ, HIV, CRP, SMQ, BMX, ALQ, DIQ tables) — no local files.

## Outputs
A single multi-page PDF, NHANES_CBC_age_trends.pdf.

## Conceptual aim / target
Characterize normal age/sex-specific reference trends for blood-count parameters in a healthy general human population.

## Note
Confirmed unrelated to the mouse CH project — a standalone human NHANES epidemiology/reference-range analysis.
