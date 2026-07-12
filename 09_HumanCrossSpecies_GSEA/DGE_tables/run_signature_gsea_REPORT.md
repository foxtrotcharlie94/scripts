# run_signature_gsea.R

**Script:** `run_signature_gsea.R`

## What it does
Builds human UP/DOWN gene "signatures" (padj<0.05) from 4 human TET2 DGE files, maps 3 mouse monocyte-subset DGE result files to human orthologs via the babelgene package, and runs fgsea testing whether each human signature is enriched in each mouse subset's ranked (ortholog-mapped) gene list, for two ranking metrics.

## Inputs
4 human DESeq2/MAST DGE files (same as run_gsea.R targets) and 3 mouse DGE CSVs (placeholder paths under DGE_tables/mouse/).

## Outputs
<mouse_subset>__<metric>.csv per subset x metric plus combined signature_GSEA_all_results.csv, written to DGE_tables/GSEA_signature_results/.

## Conceptual aim / target
Test cross-species concordance — whether the human TET2-CHIP transcriptional signature is recapitulated in specific mouse monocyte subsets' DGE rankings.

## Note
Distinct analysis — the only cross-species (mouse-to-human ortholog) signature-GSEA script in the set; shares DGE-parsing logic and human input files with run_gsea.R but answers a different question.
