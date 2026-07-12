# Shared gene-signature module (not a standalone analysis)

**Script:** `genesets_mouse.py`

## What it does
Defines the three mouse gene-symbol lists reused by nearly every other script in this project: S_PHASE and G2M_PHASE (the standard Tirosh et al. cell-cycle signature genes used by Scanpy's score_genes_cell_cycle), and PROLIFERATION_CORE (a curated ~30-gene proliferation/cell-division marker set: Mki67, Pcna, Top2a, Mcm2-7, cyclins, Cdks, Aurka/b, Plk1, etc.), bundled into a GENESETS dict.

## Inputs
None — pure data module.

## Outputs
None — imported by other scripts as `from genesets_mouse import S_PHASE, G2M_PHASE, PROLIFERATION_CORE`.

## Conceptual aim / target
Keep the cell-cycle and proliferation gene definitions consistent across every dataset/genotype analyzed in this project, so scores are directly comparable.

## Conclusions / findings
Infrastructure module — no results/conclusions of its own; underlies every proliferation and cell-cycle score reported elsewhere in this project.
