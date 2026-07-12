# =============================================================================
# quantify_coverage.R
#
# Phase-contrast cell-coverage quantification for a single image.
# Pure R, EBImage-only. Tested logic mirrors the Python demo.
#
# Usage (batch, typical case):
#   source("quantify_coverage.R")
#   quantify_coverage_folder(
#     input_dir  = "C:/data/cells/tiffs",
#     output_dir = "C:/data/cells/out"
#   )
#
# Single-file alternative:
#   quantify_coverage(
#     image_path = "C:/data/cells/well_A1.tif",
#     output_dir = "C:/data/cells/out"
#   )
#
# Outputs:
#   <output_dir>/<image_basename>/   per-image subfolder containing:
#     00_input.png            original, for reference
#     01_grayscale.png        single-channel input
#     02_flatfield.png        after white top-hat (illumination correction)
#     03_edge_density.png     Sobel gradient + Gaussian-smoothed edge map (the metric thresholded)
#     04_otsu_mask.png        raw binary mask
#     05_opened.png           after morphological opening (despeckle)
#     06_closed_filled.png    after closing + fill-hull (or copy of opened, if solidify=FALSE)
#     07_size_filtered.png    after small-object removal
#     08_overlay.png          final mask painted onto the original
#     coverage.txt            percent coverage + parameters used
#   <output_dir>/coverage_summary.csv  (batch mode only) one row per image
# =============================================================================

# --- ensure required packages are installed and loaded ----------------------
# EBImage is on Bioconductor, so we go via BiocManager (which is on CRAN).
# This block is a no-op on subsequent runs.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("EBImage", quietly = TRUE)) {
  message("Installing EBImage from Bioconductor (one-time, ~1-3 min)...")
  BiocManager::install("EBImage", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(EBImage)
})

quantify_coverage <- function(image_path,
                              output_dir,
                              edge_sigma      = 15,    # Gaussian sigma for edge-density smoothing (px)
                              bg_kernel       = 101,   # flat-field disc, must be >> cell size (odd, px)
                              open_kernel     = 3,     # opening disc (odd, px)
                              close_kernel    = 7,     # closing disc (odd, px)
                              min_object_area = 80,    # px; drop objects smaller than this
                              solidify        = TRUE,  # closing + fill-hull: fuses halo rings into solid cell bodies.
                              #   Good for sparse images; set FALSE for confluent images,
                              #   where every real gap is enclosed by cells and would be
                              #   incorrectly filled.
                              do_flatfield    = TRUE,
                              verbose         = TRUE) {
  
  stopifnot(file.exists(image_path))
  
  # ensure all brush sizes are odd (makeBrush requires odd)
  odd <- function(x) if (x %% 2L == 0L) x + 1L else as.integer(x)
  bg_kernel    <- odd(bg_kernel)
  open_kernel  <- odd(open_kernel)
  close_kernel <- odd(close_kernel)
  
  base <- tools::file_path_sans_ext(basename(image_path))
  out  <- file.path(output_dir, base)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  
  log_msg <- function(...) if (verbose) message(sprintf(...))
  
  # helper: save numeric/logical image as a viewable 8-bit PNG
  save_img <- function(im, name) {
    x <- imageData(im)
    if (is.logical(x)) x <- x * 1.0
    rng <- range(x, finite = TRUE)
    if (diff(rng) > 0) x <- (x - rng[1]) / diff(rng)
    writeImage(Image(x), file.path(out, name), type = "png")
  }
  
  # --- 0. read --------------------------------------------------------------
  log_msg("reading %s", image_path)
  raw <- readImage(image_path)
  
  # normalize() rescales to full [0,1] dynamic range for visualisation
  # (16-bit microscopy TIFFs often live in a small slice of that range)
  writeImage(normalize(raw), file.path(out, "00_input.png"), type = "png")
  
  # --- 1. grayscale ---------------------------------------------------------
  gray <- if (colorMode(raw) != Grayscale) channel(raw, "gray") else raw
  save_img(gray, "01_grayscale.png")
  
  # --- 2. flat-field (white top-hat) ---------------------------------------
  flat <- if (do_flatfield) whiteTopHat(gray, makeBrush(bg_kernel, "disc")) else gray
  save_img(flat, "02_flatfield.png")
  
  # --- 3. edge density (Sobel + Gaussian smoothing) -----------------------
  # Cells in phase-contrast generate edges from membranes, halos and internal
  # structure regardless of how confluent the field is; empty regions are
  # smooth and have near-zero edge density. This separates cleanly at any
  # confluence, unlike plain local variance which fails on confluent monolayers.
  sobel_x      <- matrix(c(-1, 0, 1,
                           -2, 0, 2,
                           -1, 0, 1), nrow = 3) / 8
  sobel_y      <- t(sobel_x)
  gx           <- filter2(flat, sobel_x)
  gy           <- filter2(flat, sobel_y)
  edges        <- sqrt(gx^2 + gy^2)
  edge_density <- gblur(edges, sigma = edge_sigma)
  save_img(edge_density, "03_edge_density.png")
  
  # --- 4. Otsu threshold on edge density -----------------------------------
  thr      <- otsu(edge_density)
  raw_mask <- edge_density > thr
  log_msg("otsu threshold = %.6f", thr)
  save_img(raw_mask, "04_otsu_mask.png")
  
  # --- 5. opening (drop salt noise) ----------------------------------------
  opened <- opening(raw_mask, makeBrush(open_kernel, "disc"))
  save_img(opened, "05_opened.png")
  
  # --- 6. closing + fill holes (solidify cells) ----------------------------
  # Skipped when solidify = FALSE so that real empty gaps between cells in a
  # confluent culture are NOT filled in.
  if (solidify) {
    closed <- closing(opened, makeBrush(close_kernel, "disc"))
    filled <- fillHull(closed)
  } else {
    filled <- opened
  }
  save_img(filled, "06_closed_filled.png")
  
  # --- 7. size filter ------------------------------------------------------
  lab <- bwlabel(filled)
  n_before <- max(lab)
  if (n_before > 0) {
    feats     <- computeFeatures.shape(lab)
    too_small <- which(feats[, "s.area"] < min_object_area)
    if (length(too_small) > 0) lab <- rmObjects(lab, too_small, reenumerate = TRUE)
  }
  final_mask <- lab > 0
  n_after    <- max(lab)
  save_img(final_mask, "07_size_filtered.png")
  
  # --- 8. overlay ----------------------------------------------------------
  overlay <- paintObjects(lab, toRGB(gray),
                          col  = c("red", "red"),
                          opac = c(1.0,  0.4))
  writeImage(overlay, file.path(out, "08_overlay.png"), type = "png")
  
  # --- coverage ------------------------------------------------------------
  coverage_pct <- mean(final_mask) * 100
  log_msg("coverage = %.3f %%  (n objects = %d, kept %d)", coverage_pct, n_before, n_after)
  
  writeLines(c(
    sprintf("image:                %s",  image_path),
    sprintf("coverage_percent:     %.4f", coverage_pct),
    sprintf("n_objects_before:     %d",  n_before),
    sprintf("n_objects_after:      %d",  n_after),
    sprintf("edge_sigma:           %g",  edge_sigma),
    sprintf("bg_kernel:            %d",  bg_kernel),
    sprintf("open_kernel:          %d",  open_kernel),
    sprintf("close_kernel:         %d",  close_kernel),
    sprintf("min_object_area:      %d",  min_object_area),
    sprintf("solidify:             %s",  solidify),
    sprintf("flatfield_correction: %s",  do_flatfield),
    sprintf("otsu_threshold:       %.6f", thr)
  ), con = file.path(out, "coverage.txt"))
  
  invisible(list(
    coverage_pct = coverage_pct,
    n_objects    = n_after,
    output_dir   = out
  ))
}

# =============================================================================
# Batch mode: process every TIFF in a folder
# =============================================================================
#
# Runs quantify_coverage() on every .tif/.tiff in input_dir. Each image gets
# its own subfolder of step PNGs + coverage.txt under output_dir, plus a
# single CSV at output_dir/<csv_name> summarising coverage_percent per image.
#
# Extra arguments (edge_sigma, solidify, etc.) are forwarded to quantify_coverage.

quantify_coverage_folder <- function(input_dir,
                                     output_dir,
                                     pattern   = "\\.tiff?$",   # case-insensitive
                                     recursive = FALSE,
                                     csv_name  = "coverage_summary.csv",
                                     ...) {
  
  stopifnot(dir.exists(input_dir))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(input_dir,
                      pattern     = pattern,
                      ignore.case = TRUE,
                      full.names  = TRUE,
                      recursive   = recursive)
  if (length(files) == 0L) stop("No TIFF files found in: ", input_dir)
  message(sprintf("Found %d TIFFs in %s", length(files), input_dir))
  
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[i]
    message(sprintf("\n[%d/%d] %s", i, length(files), basename(f)))
    res <- tryCatch(
      quantify_coverage(image_path = f, output_dir = output_dir, ...),
      error = function(e) {
        message(sprintf("  FAILED: %s", conditionMessage(e)))
        NULL
      }
    )
    rows[[i]] <- if (is.null(res)) {
      data.frame(file             = basename(f),
                 coverage_percent = NA_real_,
                 n_objects        = NA_integer_,
                 status           = "failed",
                 stringsAsFactors = FALSE)
    } else {
      data.frame(file             = basename(f),
                 coverage_percent = round(res$coverage_pct, 4),
                 n_objects        = res$n_objects,
                 status           = "ok",
                 stringsAsFactors = FALSE)
    }
  }
  
  results  <- do.call(rbind, rows)
  csv_path <- file.path(output_dir, csv_name)
  write.csv(results, csv_path, row.names = FALSE)
  message(sprintf("\nWrote summary CSV: %s   (%d / %d ok)",
                  csv_path, sum(results$status == "ok"), nrow(results)))
  
  invisible(results)
}

# -----------------------------------------------------------------------------
# Example (edit paths, then Ctrl+A, Ctrl+Enter):
# -----------------------------------------------------------------------------

res <- quantify_coverage_folder(
  input_dir  = "C:/Users/fc809/Downloads/images/mlh1cutday7/11.05.2026",
  output_dir = "C:/Users/fc809/Downloads/images/mlh1cutday7",
  solidify   = FALSE   # confluent images; set TRUE for sparse cultures
)
print(res)

# Single-file alternative:
# quantify_coverage(
#   image_path = "C:/Users/fc809/Downloads/images/mlh1ctrl10cmday7",
#   output_dir = "C:/Users/fc809/Downloads/images",
#   solidify   = FALSE
# )