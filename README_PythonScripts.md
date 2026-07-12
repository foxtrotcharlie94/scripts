# Python Scripts Portfolio — Index

Organized by project. Each script has a matching `<name>_REPORT.md` / `<name>_REPORT.docx` explaining what it does, its inputs/outputs, its conceptual aim, and its conclusions.


## 01_BMDM_Myeloid_Inflammation

3 scripts. See `_PROJECT_OVERVIEW.md` in that folder for the project-level summary.

- `bmdm_asymmetry_diagnostics.py` — BMDM library-composition & logFC tail-asymmetry diagnostics
- `bmdm_axis_concordance.py` — Burden-axis vs. lenti-axis concordance analysis
- `ch_myeloid_heatmap.py` — CH / myeloid-inflammation gene panel heatmap

## 02_Division_Burden_Methylation_Clocks

8 scripts. See `_PROJECT_OVERVIEW.md` in that folder for the project-level summary.

- `analyze_gse44117.py` — GSE44117 global methylation by group (Fig 2E setup)
- `apply_clocks.py` — Apply Blood + WLMT epigenetic clocks to GSE44117 RRBS data
- `apply_clocks2.py` — Apply YOMT + Liver clocks and compute per-CpG age/reconstitution deltas
- `cpg_delta_correlation.py` — Core hypothesis test: aging delta vs. reconstitution delta per CpG (Blood, WLMT)
- `cpg_delta_overlay.py` — Overlay figure: Blood + WLMT correlation on one panel (v1)
- `cpg_delta_overlay3.py` — Overlay figure: Blood + WLMT + YOMT correlation (v2, 3 clocks)
- `cpg_delta_overlay4.py` — Overlay figure: Blood + WLMT + YOMT correlation (v3, recolored)
- `plot_clock_scores.py` — Weighted clock scores by group, plus positive/negative-weight CpG split

## 03_Mutant_HSPC_Proliferation_MultiGenotype

29 scripts. See `_PROJECT_OVERVIEW.md` in that folder for the project-level summary.

- `genesets_mouse.py` — Shared gene-signature module (not a standalone analysis)
- `bulk_analysis.py` — First-pass bulk RNA-seq proliferation signature scoring (3 CH drivers)
- `rigor_bulk.py` — Rigorous re-analysis: DESeq2 + competitive proliferation gene-set test (all 3 datasets)
- `rigor_one.py` — Rigorous DESeq2 test for a single dataset (Asxl1 or Jak2)
- `s1_smartseq2.py` — Smart-seq2 LT-HSC cell-cycle scoring (WT vs Tet2 vs Dnmt3a)
- `s2_tenx.py` — 10x LSK cell-cycle scoring, single genotype (GSE263101, WT or Jak2V617F)
- `s3_synth.py` — Cross-dataset synthesis: phase fractions & proliferation across all 3 mutants
- `s4_clusters_ss2.py` — Smart-seq2: add Leiden clustering + lineage annotation to WT/Tet2/Dnmt3a
- `scrna_analysis.py` — Consolidated single-cell cell-cycle analysis (Smart-seq2 + 10x, one script)
- `c1_build_ss2.py` — Smart-seq2 pipeline step 1: build merged, scored, PCA'd dataset
- `c2_cluster_ss2.py` — Smart-seq2 pipeline step 2: Leiden clustering + per-cluster proliferation
- `c_ss2_full.py` — Smart-seq2 build+cluster in one script (KMeans fallback)
- `c_tenx_celltype.py` — 10x single-genotype cell-type + proliferation scoring (Jak2 dataset)
- `c_join.py` — 10x WT vs Jak2V617F: compartment-matched proliferation comparison
- `c_final.py` — Final cross-genotype summary figure (matched compartment)
- `audit_ss2.py` — QC audit: does sequencing depth confound the Smart-seq2 proliferation signal?
- `ax_panel.py` — Asxl1 10x panel-extraction step (per sample tag)
- `ax_join.py` — Asxl1 10x: joined compartment-matched WT vs KO comparison (wk04 + wk36)
- `asxl1_10x.py` — Asxl1 10x: self-contained single-script version of ax_panel.py + ax_join.py
- `axd_panel.py` — Asxl1 10x panel extraction with in-matrix hashtag demultiplexing
- `d_panel.py` — Dnmt3a 10x panel-extraction step (per sample, GSE272266)
- `jak2_panel.py` — Jak2 (GSE227026) 10x panel-extraction step (per sample)
- `r878_panel.py` — Dnmt3a-R878H 10x panel-extraction step with full-lineage marker panel
- `t_panel.py` — Jak2 (GSE263101) 10x panel-extraction step used by c_join.py
- `t_save.py` — Jak2 (GSE263101) 10x: save full log-normalized matrix (all genes)
- `t2_panel.py` — Tet2 (GSE209994) 10x panel-extraction step with HTO demultiplexing
- `t2h_panel.py` — Tet2 (GSE209994) 10x panel-extraction step — near-duplicate of t2_panel.py
- `strict_join.py` — Generalized, stricter compartment-matched join/comparison script
- `g2m_vs_s_signtest.py` — Asxl1 reconciliation: is the proliferation signal S-phase or G2M-driven?
