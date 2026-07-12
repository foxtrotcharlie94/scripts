# lenti_overlaps_with_mouse_genome.R

**Script:** `lenti_overlaps_with_mouse_genome.R`

## What it does
BLASTs a lentiviral vector FASTA against a custom-built mouse mm10 transcriptome database (mRNA, lncRNA, pre-mRNA/gene-body sequences), filters hits by identity/length/e-value, annotates hit gene symbol and transcript region (CDS/UTR/intron), and produces summary plots and masked/deleted versions of the vector sequence.

## Inputs
A vector FASTA (my_vector.fasta) and mm10 genome/annotation packages used to build the BLAST DB.

## Outputs
vector_transcriptome_overlaps_full.csv, vector_transcriptome_overlaps_summary.csv, overlap_summary_plot.pdf, vector_overlaps_deleted.fasta, vector_overlaps_Nmasked.fasta.

## Conceptual aim / target
Identify which parts of the lentiviral clonal-marking vector share homology with mouse transcripts, to flag off-target mapping/cross-hybridization risk.

## Note
Distinct analysis — unrelated in mechanism to the other scripts, though motivated by the same lentiviral clonal-marking project.
