# First-pass bulk RNA-seq proliferation signature scoring (3 CH drivers)

**Script:** `bulk_analysis.py`

**Status:** First-pass / superseded

## What it does
Loads three independent public bulk RNA-seq datasets for Tet2 (GSE96758), Asxl1 (GSE158184) and Jak2V617F (GSE123401), z-scores log-CPM expression per gene, scores each sample against the S-phase/G2M/proliferation-core gene sets, and runs a per-gene-set Welch t-test comparing mutant vs WT samples in each dataset.

## Inputs
TET2GSE96758_normalizedData.txt.gz, GSE158184_Bulk_ASXL_counts.filt.txt.gz, and GSE123401_RAW.tar + GSE263101_RAW.tar (the latter supplies an Ensembl-to-symbol map for the featureCounts-format Jak2 data).

## Outputs
<dataset>_signature_scores_persample.csv (per sample, per gene set) and bulk_proliferation_summary.csv (t-test results per dataset per gene set).

## Conceptual aim / target
Quick first-pass check across three unrelated CH driver mutations of whether cell-cycle/proliferation gene-set expression is shifted in mutant vs WT bulk samples, before moving to more rigorous statistics.

## Conclusions / findings
Superseded by rigor_bulk.py / rigor_one.py, which replace the ad-hoc z-score/t-test approach with proper DESeq2 differential expression plus competitive gene-set testing (this naive version doesn't correct for the fact that gene-set means are not independent observations).
