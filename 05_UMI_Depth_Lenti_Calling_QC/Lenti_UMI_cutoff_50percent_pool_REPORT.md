# Lenti_UMI_cutoff_50percent_pool.R

**Script:** `Lenti_UMI_cutoff_50percent_pool.R`

## What it does
Computes CDS (coding sequence) sizes and merged exon coordinates for two specific mouse genes (Kmt2d and Pkd1) using TxDb.Mmusculus.UCSC.mm10.knownGene annotation — maps gene symbols to Entrez IDs, retrieves per-transcript CDS intervals, merges overlapping exon blocks, and totals CDS length per gene.

## Inputs
Bioconductor annotation packages TxDb.Mmusculus.UCSC.mm10.knownGene and org.Mm.eg.db (no user data files).

## Outputs
cds_sizes_Kmt2d_Pkd1.csv and cds_coordinates_Kmt2d_Pkd1.bed.

## Conceptual aim / target
Obtain accurate CDS length/coordinates for Kmt2d and Pkd1, likely to normalize mutation/variant calling or coverage by gene size.

## Note
Filename is misleading — despite suggesting a Lenti-UMI-cutoff analysis, the actual code is an unrelated CDS-size/BED-coordinate extraction script for Kmt2d and Pkd1.
