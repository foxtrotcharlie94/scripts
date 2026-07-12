# Project overview — H2B Label-Retention, Physiological Aging & CytoTRACE2

Tests whether the lab's HB/LB and LentiPos/Neg LT-HSC/monocyte signatures relate to two external
reference axes: (1) an H2B-GFP label-retention (quiescent vs dividing HSC) signature
(`H2B_GSEA/H2B_analysis.R`, `H2B_fdrcutoff/H2B_analysis.R` — a fixed-N vs FDR-based gene-set-
definition sensitivity pair), and (2) physiological monocyte aging from the public GSE207063
dataset (`AgingCell_analysis.R`, and its much larger, repeatedly-iterated sibling
`H2B_fdrcutoff/H2B_pathways_wFDR.R`, which also folds in a final TET2-signature comparison).
`cytotrace2_LT-HSCs.R` uses the third-party CytoTRACE2 package plus Slingshot pseudotime to test
whether burden/marking status shifts LT-HSCs toward a more differentiated state.
