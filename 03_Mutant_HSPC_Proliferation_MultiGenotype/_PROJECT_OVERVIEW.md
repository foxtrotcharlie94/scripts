# Project overview — Mutant HSPC Proliferation, R pipeline (Seurat/Harmony)

These four R scripts are the single-cell clustering backbone behind the Python proliferation
comparison project: `run_seurat_harmony_all_mutants.R` builds a harmonized, clustered HSPC atlas
across five CH driver genotypes (Tet2, Asxl1, Jak2, Dnmt3a, Dnmt3a-R878H) from public 10x
datasets; `cluster_markers.R` annotates the resulting clusters with lineage markers;
`fix_cluster_markers.R` and `finish_asxl1_res1.2_markers.R` are bugfix/patch reruns after the
original marker-finding step silently dropped large progenitor clusters. The compartment-matched
WT-vs-mutant proliferation comparisons themselves (the `c_join.py`/`c_final.py` results) are
Python, downstream of these R-built Seurat objects.
