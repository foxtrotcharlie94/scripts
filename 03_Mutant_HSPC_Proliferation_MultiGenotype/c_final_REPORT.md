# Final cross-genotype summary figure (matched compartment)

**Script:** `c_final.py`

## What it does
Hardcodes the final, matched-compartment result for all three genotypes analyzed with single-cell data (Tet2 LT-HSC Smart-seq2: WT 42.4% vs MUT 28.8% cycling, p=0.0076, down; Dnmt3a LT-HSC Smart-seq2: 42.4% vs 45.3%, p=0.74, not significant; Jak2V617F LT-HSC+MPP1 10x: 45.7% vs 52.9%, p=5.5e-19, up) and draws the project's headline summary bar chart.

## Inputs
None loaded from disk — the three result rows are hardcoded from the upstream cluster/compartment analyses (c2_cluster_ss2.py / c_ss2_full.py for Tet2 & Dnmt3a, c_join.py for Jak2V617F).

## Outputs
MASTER_matched_compartment.csv, fig6_summary.png.

## Conceptual aim / target
Aggregate the project's three genotype-specific, compartment-matched proliferation results into one directly comparable summary figure and table.

## Conclusions / findings
The project's headline conclusion: matching genotype to its correct primitive compartment reveals genotype-specific, even opposite-direction effects on proliferation — Tet2-KO LT-HSCs cycle LESS than WT, Dnmt3a-mutant LT-HSCs show no significant difference, and Jak2V617F LT-HSC+MPP1 cells cycle substantially MORE than WT — underscoring that pooled, non-compartment-matched comparisons (as in the earlier s1-s4/scrna_analysis.py scripts) would have obscured or reversed these effects.
