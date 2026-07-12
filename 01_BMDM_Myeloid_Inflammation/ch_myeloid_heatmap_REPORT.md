# CH / myeloid-inflammation gene panel heatmap

**Script:** `ch_myeloid_heatmap.py`

## What it does
Builds a two-panel heatmap over an 18-gene curated panel spanning cytokines (Il1b, Il6, Tnf...), inflammasome components (Nlrp3, Casp1...), chemokines (Ccl2, Cxcl1...) and signaling/effector genes (Tlr2, Nfkb1, Nos2...). Panel A shows row-z-scored mean log2 CPM expression per group; Panel B shows the signed effect size (log2 mean-CPM difference) for four pairwise contrasts with BH-FDR significance stars — burden contrasts use an unpaired Welch t-test, lenti contrasts use a donor-paired t-test.

## Inputs
Same DEG/CPM files as the other two BMDM scripts; a hardcoded 18-gene panel list is defined in the script.

## Outputs
ch_myeloid_heatmap.png, ch_myeloid_heatmap_stats.csv (mean CPM per group, log2 difference, p-value and BH-q per contrast per gene); console list of genes reaching BH-q<0.05 per contrast.

## Conceptual aim / target
Directly visualize, with proper statistical pairing, whether specific inflammatory/CH-associated genes are induced by burden (mutant clone size), by the lentiviral construct itself, or both.

## Conclusions / findings
Identifies which individual inflammatory genes reach BH-q<0.05 in each of the four contrasts (printed at runtime and saved in the stats CSV); used to pinpoint specific candidate genes rather than only genome-wide trends.
