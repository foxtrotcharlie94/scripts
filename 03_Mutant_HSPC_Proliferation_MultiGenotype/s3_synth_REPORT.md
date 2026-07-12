# Cross-dataset synthesis: phase fractions & proliferation across all 3 mutants

**Script:** `s3_synth.py`

**Status:** Superseded by c_final.py (compartment-matched version)

## What it does
Combines the per-cell outputs of s1_smartseq2.py (WT/Tet2/Dnmt3a) and s2_tenx.py (WT/Jak2V617F) into master summary tables and three comparison figures: stacked cell-cycle-phase composition bars, a cycling-fraction bar chart, and per-cell proliferation-score violin plots, all spanning the three CH driver genotypes and their respective WT controls.

## Inputs
obs_smartseq2_WT/Tet2/Dnmt3a.csv and obs_tenx_WT/Jak2V617F.csv (outputs of s1_smartseq2.py and s2_tenx.py).

## Outputs
MASTER_scrna_phase_fractions.csv, MASTER_scrna_comparisons.csv, fig1_phase_fractions.png, fig2_cycling_fraction.png, fig3_prolif_score.png.

## Conceptual aim / target
Produce one unified, cross-genotype summary of the pooled (non-compartment-matched) single-cell proliferation comparison, as an interim checkpoint before adding cluster/compartment structure.

## Conclusions / findings
The MASTER_scrna_comparisons.csv table is the pooled-cell-level precursor to the more rigorous, compartment-matched comparison later finalized in c_final.py's MASTER_matched_compartment.csv.
