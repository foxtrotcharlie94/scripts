# Dnmt3a 10x panel-extraction step (per sample, GSE272266)

**Script:** `d_panel.py`

**Status:** Panel-extraction step; no join/comparison script recovered

## What it does
Per-sample preprocessing step for the Dnmt3a GSE272266 10x dataset (replicated across vehicle/sex/hormone arms). Extracts one sample's matrix/features/barcodes from the dataset's tar archive by GEO sample ID, QC-filters, normalizes, scores cell cycle + proliferation-core, and either defines or reuses a shared gene panel — the Dnmt3a analogue of ax_panel.py.

## Inputs
GSE272266_RAW.tar; GEO sample ID, genotype, sex, replicate number and mode passed as command-line arguments.

## Outputs
/tmp/dnmt3a/panel.csv (shared panel) and /tmp/dnmt3a/<geno>_<sex><rep>_X.npy/_obs.csv per sample.

## Conceptual aim / target
Prepare per-sample, per-cell scored expression panels for the more highly replicated Dnmt3a dataset identified in BETTER_DATASETS_Tet2_Dnmt3a.md, following the same panel-extraction template used for Asxl1, as a step toward a compartment-matched Dnmt3a WT-vs-mutant comparison with sex as a covariate.

## Conclusions / findings
Preparation step only; per-sample cell counts printed at runtime. No corresponding "d_join.py" was found among the recovered scripts — this upgraded-dataset Dnmt3a comparison appears to have been left at the panel-extraction stage (the Dnmt3a conclusion ultimately reported in c_final.py instead uses the earlier, less-replicated Smart-seq2 GSE124822 data).
