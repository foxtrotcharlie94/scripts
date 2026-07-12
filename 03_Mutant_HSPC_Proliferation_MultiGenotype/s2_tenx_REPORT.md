# 10x LSK cell-cycle scoring, single genotype (GSE263101, WT or Jak2V617F)

**Script:** `s2_tenx.py`

**Status:** Superseded by c_tenx_celltype.py / c_join.py

## What it does
10x-data counterpart to s1_smartseq2.py: loads one genotype's sample (WT or Jak2V617F) from the GSE263101 10x tar archive, filters low-quality/high-mitochondrial cells, normalizes, and scores cell-cycle phase + proliferation-core signature per cell. Genotype and file prefix are passed as command-line arguments so it is run once per genotype.

## Inputs
GSE263101_RAW.tar (10x matrix/features/barcodes for a given sample prefix).

## Outputs
obs_tenx_<genotype>.csv (per-cell phase/scores); console summary of phase percentages and mean proliferation score.

## Conceptual aim / target
Generate the per-cell cell-cycle/proliferation scores needed to compare WT vs Jak2V617F LSK cells, matching the same scoring methodology used for the Smart-seq2 genotypes.

## Conclusions / findings
Feeds directly into s3_synth.py's cross-dataset comparison; superseded for compartment-level conclusions by the c_join.py / c_tenx_celltype.py pipeline.
