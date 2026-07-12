# Well_coverage_drightfield_analysis.R

**Script:** `Well_coverage_drightfield_analysis.R`

## What it does
Pure-R (EBImage-based) phase-contrast cell image analysis pipeline: reads a TIFF, converts to grayscale, applies flat-field correction, computes an edge-density map (Sobel + Gaussian smoothing), Otsu-thresholds it, cleans up morphologically, then quantifies percent cell coverage; includes a batch-mode wrapper over a folder of images.

## Inputs
TIFF/TIF brightfield/phase-contrast well images.

## Outputs
Per-image step-by-step PNGs (grayscale, flatfield, edge density, masks, overlay), coverage.txt, plus batch coverage_summary.csv.

## Conceptual aim / target
Quantify percent confluence/cell coverage in brightfield microscopy wells — a cell-culture imaging QC tool.

## Note
Unrelated to the CH/lentiviral single-cell project — a brightfield image-analysis (confluence quantification) script.
