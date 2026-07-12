# 10x single-genotype cell-type + proliferation scoring (Jak2 dataset)

**Script:** `c_tenx_celltype.py`

**Status:** Alternate (marker-argmax) cell-typing approach, feeds c_join.py

## What it does
Extracts one genotype's sample from the GSE263101 10x tar archive, QC-filters, normalizes, scores cell cycle and proliferation-core signature, then scores 7 lineage-marker sets (HSC, MPP, MkP, Ery, GMP/Mye, Lymph, Baso/MC) per cell and assigns each cell to its single top-scoring lineage — a per-cell cell-typing step rather than cluster-based typing.

## Inputs
GSE263101_RAW.tar; genotype label and file prefix passed as command-line arguments.

## Outputs
tenx_celltype_<genotype>.csv (per-cell genotype, assigned celltype, phase, cycling flag, proliferation score); console cell-type counts.

## Conceptual aim / target
Provide a lineage label for every 10x cell using direct marker-score argmax (rather than clustering first), as a simpler/faster alternative cell-typing approach for the Jak2 10x dataset.

## Conclusions / findings
Produces per-cell lineage calls used as one input path toward the matched-compartment comparison; the cluster-based approach in c_join.py was used for the actual statistical comparison reported in c_final.py.
