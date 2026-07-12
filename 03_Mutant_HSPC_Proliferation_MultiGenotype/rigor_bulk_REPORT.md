# Rigorous re-analysis: DESeq2 + competitive proliferation gene-set test (all 3 datasets)

**Script:** `rigor_bulk.py`

## What it does
Redoes the bulk proliferation comparison properly: runs DESeq2 for Asxl1 (4v4) and Jak2V617F (3v3), and a Welch t-test on pre-normalized data for Tet2 (3v3, no raw counts available). For each dataset, tests whether proliferation-gene log fold-changes are competitively shifted relative to the rest of the genome (Mann-Whitney U) and whether they are significantly different from zero as a group (Wilcoxon signed-rank), plus reports key single-gene logFC/padj values (Mki67, Top2a, etc.).

## Inputs
GSE158184_Bulk_ASXL_counts.filt.txt.gz, GSE263101_RAW.tar + GSE123401_RAW.tar (Jak2), TET2GSE96758_normalizedData.txt.gz.

## Outputs
bulk_proliferation_RIGOROUS.csv (one row per dataset: median LFC proliferation vs genome, competitive p, signed Wilcoxon p, fraction of proliferation genes up).

## Conceptual aim / target
Replace the naive per-gene-set z-score t-test in bulk_analysis.py with a statistically defensible test that (a) uses real differential-expression modeling (DESeq2) where raw counts exist, and (b) tests proliferation genes competitively against genome-wide background rather than in isolation.

## Conclusions / findings
Produces the per-dataset competitive-test p-values and fraction of proliferation genes up/down that feed into the project's overall cross-genotype comparison (see MASTER_scrna_comparisons.csv / c_final.py for how this integrates with the single-cell results).
