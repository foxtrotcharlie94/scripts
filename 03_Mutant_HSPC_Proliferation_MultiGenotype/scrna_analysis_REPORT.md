# Consolidated single-cell cell-cycle analysis (Smart-seq2 + 10x, one script)

**Script:** `scrna_analysis.py`

**Status:** Alternate all-in-one version of s1+s2+s3

## What it does
A self-contained alternative that combines what s1_smartseq2.py + s2_tenx.py + s3_synth.py do separately: processes the Smart-seq2 (WT/Tet2/Dnmt3a) and 10x (WT/Jak2V617F) datasets in one run, scores cell cycle/proliferation for each, and computes the same chi-square/Mann-Whitney WT-vs-mutant comparisons for all three genotypes without writing intermediate per-genotype CSVs first.

## Inputs
Same GSE124822 Smart-seq2 files and GSE263101_RAW.tar as s1/s2.

## Outputs
scrna_cellcycle_summary.csv, scrna_comparisons.csv.

## Conceptual aim / target
Provide a single, reproducible entry point for the pooled-cell cross-genotype proliferation comparison, functionally equivalent to running s1+s2+s3 in sequence.

## Conclusions / findings
Produces the same conclusions as the s1/s2/s3 chain (pooled-cell WT-vs-mutant cycling/proliferation differences per genotype); superseded for compartment-matched analysis by the c_*/strict_join.py pipeline.
