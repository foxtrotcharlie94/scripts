# QC audit: does sequencing depth confound the Smart-seq2 proliferation signal?

**Script:** `audit_ss2.py`

## What it does
Re-processes the same Smart-seq2 plates (WT/Tet2_30/Tet2_272/Dnmt3a_675, plus a second WT plate) with per-plate/batch QC summaries (median genes detected, median total counts, %cycling), checks the within-WT correlation between genes-detected and cell-cycle score, and fits depth-adjusted logistic regression models (cycling ~ genotype + sequencing-depth z-score) for Tet2 and Dnmt3a vs WT to see whether the naive (unadjusted) proliferation difference survives controlling for detection depth.

## Inputs
Same GSE124822 Smart-seq2 RNA count files as s1_smartseq2.py/c1_build_ss2.py.

## Outputs
ss2_audit_percell.csv, ss2_audit_perplate.csv; console odds ratios (unadjusted vs. depth-adjusted) and p-values per genotype.

## Conceptual aim / target
Rule out that any apparent Tet2/Dnmt3a proliferation shift is actually an artifact of different sequencing depth between plates/batches (deeper sequencing tends to detect more transcripts, inflating cell-cycle scores).

## Conclusions / findings
Reports, per mutant, both the unadjusted and depth-adjusted odds ratio and p-value for cycling ~ genotype — used to confirm whether the c2/c_ss2_full cluster-level proliferation differences hold up after controlling for detection depth as a technical confound.
