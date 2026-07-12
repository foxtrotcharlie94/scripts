# combine_barcodes.R
# Combine per-sample barcode/UMI count files into a single wide matrix.
# Each input file is a 2-column TSV (no header): barcode, UMI_count.
# Output: one row per unique barcode, one column per sample, UMI counts in cells.
# Barcodes absent from a sample are filled with 0.

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- USER SETTINGS ------------------------------------------------------
# Folder containing all the per-sample *_Barcode_UMI_number.txt files
input_dir   <- "."                              # change to your folder
output_file <- "combined_barcodes_wide.tsv"

# Strip the boilerplate suffix to derive a sample name from the filename.
# e.g. "5_Blood_Granulocytes_CA_S82_Barcode_UMI_number.txt" -> "5_Blood_Granulocytes_CA_S82"
filename_suffix <- "_Barcode_UMI_number\\.txt$"

# Drop the "[]" pseudo-entry (unmapped / no-barcode reads)?
drop_empty_barcode <- TRUE
# -------------------------------------------------------------------------

files <- list.files(
  input_dir,
  pattern     = filename_suffix,
  full.names  = TRUE
)

if (length(files) == 0) {
  stop("No files matched pattern '", filename_suffix, "' in '", input_dir, "'")
}

message("Found ", length(files), " files.")

read_barcode_file <- function(path) {
  sample_name <- sub(filename_suffix, "", basename(path))
  readr::read_tsv(
    path,
    col_names      = c("barcode", "umi_count"),
    col_types      = "ci",
    show_col_types = FALSE
  ) |>
    dplyr::mutate(sample = sample_name)
}

long <- purrr::map_dfr(files, read_barcode_file)

# Defensive: if any (sample, barcode) appears more than once, sum the UMIs.
long <- long |>
  dplyr::group_by(sample, barcode) |>
  dplyr::summarise(umi_count = sum(umi_count), .groups = "drop")

if (drop_empty_barcode) {
  long <- dplyr::filter(long, barcode != "[]")
}

# Pivot to wide: rows = barcodes, cols = samples, values = UMI counts (0 if absent)
wide <- long |>
  tidyr::pivot_wider(
    names_from  = sample,
    values_from = umi_count,
    values_fill = 0
  )

# Add a couple of convenience columns: number of samples each barcode is detected in,
# and total UMIs across all samples. Move barcode + summary cols to the front.
sample_cols <- setdiff(colnames(wide), "barcode")
wide <- wide |>
  dplyr::mutate(
    n_samples_detected = rowSums(dplyr::across(dplyr::all_of(sample_cols)) > 0),
    total_umi          = rowSums(dplyr::across(dplyr::all_of(sample_cols)))
  ) |>
  dplyr::relocate(barcode, n_samples_detected, total_umi)

# Sort by total abundance (descending), then by barcode for ties
wide <- dplyr::arrange(wide, dplyr::desc(total_umi), barcode)

readr::write_tsv(wide, output_file)

message("Wrote ", output_file)
message("Unique barcodes: ", nrow(wide))
message("Samples:         ", length(sample_cols))
