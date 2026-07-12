# neutral_region_pipeline.R

**Script:** `neutral_region_pipeline.R`

## What it does
Identifies genomic "neutral" candidate regions in mm10 for a liver mutagenesis assay by building an exclusion mask (ENCODE blacklist, gene bodies+buffer, cCREs, repeats, segmental duplications, imprinted regions, liver ChIP-seq peaks), tiling into 200kb windows, scoring by GC content/CpG O:E/mappability, then ranking/selecting top windows per chromosome; a second section chops top candidates into 225bp fragments for assay design.

## Inputs
mm10-blacklist.v2.bed.gz, mm10-cCREs.bed.gz, rmsk.txt.gz, genomicSuperDups.txt.gz, mm10.k100.umap.bedgraph.gz, mm10_imprinted_regions.bed, several liver ChIP-seq peak BED files.

## Outputs
exclusion_mask_mm10_liver.bed, candidates_ranked_full.tsv, top_candidates_liver.bed, top_candidates_usable_fragments.bed, assay_fragments_225bp.bed/_full.tsv.

## Conceptual aim / target
Find neutral, non-functional, mappable genomic loci in mouse liver suitable as safe-harbor/negative-control regions for a mutagenesis/barcoding assay.

## Note
Distinct analysis — unrelated to the CH/GSEA scripts; a liver-tissue genome-annotation/target-selection pipeline.
