# Project overview — Multi-Genotype HSPC Proliferation Comparison

**Biological question:** across several unrelated clonal-hematopoiesis (CH)
driver mutations (Asxl1-KO, Tet2-KO, Dnmt3a-mutant, Jak2-V617F), is
hematopoietic stem/progenitor cell (HSPC) proliferation increased, and is any
such effect specific to the primitive stem/progenitor compartment rather than
more differentiated cells?

**Design:** this is the largest and most iterative project in the portfolio —
28 scripts spanning bulk RNA-seq (3 independent public datasets), Smart-seq2
single-cell (GSE124822: WT/Tet2/Dnmt3a), and 10x single-cell (GSE263101 Jak2,
GSE158184 Asxl1, GSE209994 Tet2, GSE272266 Dnmt3a, plus exploratory GSE227026/
R878H datasets), scored throughout with a shared Tirosh cell-cycle +
proliferation-core gene signature (`genesets_mouse.py`).

**How to read this folder:** the scripts fall into four rough generations —

1. **First-pass / pooled-cell scoring** (`bulk_analysis.py`, `s1_smartseq2.py`,
   `s2_tenx.py`, `s3_synth.py`, `s4_clusters_ss2.py`, `scrna_analysis.py`) —
   naive z-score/chi-square comparisons pooling all cells per genotype.
2. **Rigorous re-analysis** (`rigor_bulk.py`, `rigor_one.py`, `audit_ss2.py`)
   — proper DESeq2 + competitive gene-set tests, and a depth-confound audit.
3. **Compartment-matched pipeline** (`c1_build_ss2.py` → `c2_cluster_ss2.py`,
   `c_ss2_full.py`, `c_tenx_celltype.py`, `c_join.py`, `strict_join.py`) —
   clusters cells first and restricts the WT-vs-mutant comparison to the
   correctly matched primitive (LT-HSC+MPP1) compartment.
4. **Per-genotype panel-extraction scripts** (`ax_panel.py`/`ax_join.py`/
   `asxl1_10x.py`/`axd_panel.py` for Asxl1; `d_panel.py` for Dnmt3a;
   `jak2_panel.py`/`t_panel.py`/`t_save.py` for Jak2; `t2_panel.py`/
   `t2h_panel.py` for Tet2; `r878_panel.py` for Dnmt3a-R878H) — a reusable
   template for loading one genotype's 10x sample, scoring it, and extracting
   a shared gene panel for later joining.

**Headline result:** `c_final.py` aggregates the three genotypes that reached
a compartment-matched conclusion — Tet2-KO LT-HSCs cycle significantly LESS
than WT (p=0.0076), Dnmt3a-mutant LT-HSCs show no significant difference
(p=0.74), and Jak2V617F LT-HSC+MPP1 cells cycle significantly MORE than WT
(p=5.5e-19) — demonstrating genotype-specific, non-uniform effects on HSPC
proliferation that would have been missed or reversed by pooled,
non-compartment-matched analysis.

`g2m_vs_s_signtest.py` (in the `seurat_harmony_output/` subfolder) is a later
reconciliation step against a separate R/Seurat clustering of Asxl1 data,
showing that an apparent Asxl1-KO cell-cycle increase was not a clean,
symmetric S+G2M proliferation increase once the two phases were examined
separately — a caution about over-interpreting combined cell-cycle scores.

Several panel-extraction scripts (`d_panel.py`, `jak2_panel.py`,
`r878_panel.py`) have no corresponding "join"/comparison script recovered —
these upgraded/exploratory datasets appear to have been prepared but not
carried through to a final statistical comparison.
