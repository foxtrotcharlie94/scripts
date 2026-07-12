# Smart-seq2 LT-HSC cell-cycle scoring (WT vs Tet2 vs Dnmt3a)

**Script:** `s1_smartseq2.py`

**Status:** Superseded by c1_build_ss2.py / c2_cluster_ss2.py

## What it does
Loads the GSE124822 Smart-seq2 LT-HSC plate count data for WT, Tet2-KO and Dnmt3a-mutant genotypes, normalizes/log-transforms, scores each cell for Tirosh cell-cycle phase and the proliferation-core signature, summarizes %cells in each phase per genotype, and tests WT vs each mutant (chi-square on cycling-vs-not counts; Mann-Whitney U on proliferation score).

## Inputs
GSE124822_LTHSC-WT_476/846, GSE124822_LTHSC-TET2_30/272, GSE124822_LTHSC-DNMT3A_675 RNA count text files.

## Outputs
scrna_smartseq2_summary.csv (phase fractions + proliferation mean per genotype) and obs_smartseq2_<genotype>.csv (per-cell scores); console chi-square/MWU results.

## Conceptual aim / target
First single-cell-resolution look at whether Tet2-KO and Dnmt3a-mutant LT-HSCs show a different proliferating fraction than WT, pooling all profiled LT-HSCs per genotype.

## Conclusions / findings
Provides the pooled (non-compartment-matched) WT-vs-mutant cycling % and proliferation-score comparison later superseded by the compartment-matched c1/c2/c_ss2_full pipeline, which separates true LT-HSC-like cells from more differentiated cells within the same plates.
