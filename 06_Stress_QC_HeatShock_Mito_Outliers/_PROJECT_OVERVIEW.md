# Project overview — Stress/Heat-Shock/Mitochondrial QC & Outlier Characterization

Scripts investigating whether particular samples (especially HB3, HB4, HB2, HB5) are technical
or biological outliers, and whether their divergent transcriptomes reflect a cellular stress
response (heat-shock, Trp53, mitochondrial/OXPHOS) rather than the biology under study. Includes
both a "raw fraction" and a "z-scored module score" approach to quantifying per-cell heat-shock
signal (for both LT-HSCs and classical monocytes), plus dedicated outlier-sample PCA/volcano/ORA
analyses (`HB34_v_HB125_analysis.R`, `H2_H5_analysis.R`).
