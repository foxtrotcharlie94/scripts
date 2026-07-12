# barcode_visualizations.R
# Reads combined_barcodes_wide.tsv and produces:
#   1.  Histogram of barcode sharing across samples
#   2.  Barcodes per sample, faceted by mouse
#   3a. UpSet plot per mouse (all samples)
#   3b. Venn diagram per mouse PER hematopoietic POPULATION
#       (sets within each Venn = the different bone/source samples)
#   4.  Heatmap of barcode counts (mouse x sample_type)
# Each plot saved as PNG and PDF; summary tables saved as CSV.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)   # install.packages("ggVennDiagram")
  library(ComplexUpset)    # install.packages("ComplexUpset")
  library(scales)
})

# ---- USER SETTINGS ------------------------------------------------------
input_file <- "combined_barcodes_wide.tsv"   # in working dir, or give full path
output_dir <- "C:/Users/fc809/Downloads/C_regions_results/"

# Hematopoietic populations to make per-mouse Venns for. The script will draw
# one Venn per (mouse x population) combination, with each set = one bone/source
# sample for that mouse+population. Combos with <2 samples are skipped.
populations <- c("GMP", "MPP4", "LSK", "CLP")
# -------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# ---- LOAD & PARSE METADATA ---------------------------------------------
wide <- readr::read_tsv(input_file, show_col_types = FALSE)

summary_cols <- c("barcode", "n_samples_detected", "total_umi")
sample_cols  <- setdiff(colnames(wide), summary_cols)

# Parse sample columns of the form: {mouse}_{population}_{source}_S{number}
# e.g. 5_LSK_Spine_CA_S62 -> mouse=5, population=LSK, source=Spine_CA, sample_type=LSK_Spine_CA
meta <- tibble::tibble(sample = sample_cols) |>
  dplyr::mutate(
    mouse       = as.integer(stringr::str_extract(sample, "^\\d+")),
    population  = stringr::str_extract(sample, "(?<=^\\d{1,3}_)[^_]+"),
    source      = sample |>
      stringr::str_remove("^\\d+_[^_]+_") |>   # strip "{mouse}_{population}_"
      stringr::str_remove("_S\\d+$"),           # strip trailing "_S##"
    sample_type = sample |>
      stringr::str_remove("^\\d+_") |>
      stringr::str_remove("_S\\d+$")
  )

if (any(is.na(meta$mouse))) {
  warning("Some sample columns did not start with a digit; check column naming.")
}
message("Parsed ", nrow(meta), " samples across ",
        dplyr::n_distinct(meta$mouse), " mice (",
        paste(sort(unique(meta$mouse)), collapse = ", "), ").")
message("Populations detected: ",
        paste(sort(unique(meta$population)), collapse = ", "))

# Long form for downstream tallies
long <- wide |>
  dplyr::select(barcode, all_of(sample_cols)) |>
  tidyr::pivot_longer(-barcode, names_to = "sample", values_to = "umi_count") |>
  dplyr::filter(umi_count > 0) |>
  dplyr::left_join(meta, by = "sample")

sample_type_levels <- sort(unique(meta$sample_type))

# ---- 1. HISTOGRAM -------------------------------------------------------
p1 <- wide |>
  dplyr::count(n_samples_detected) |>
  ggplot(aes(n_samples_detected, n)) +
  geom_col(fill = "#4C72B0") +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  labs(x = "Number of samples a barcode is detected in",
       y = "Number of barcodes",
       title = "Barcode sharing across samples") +
  theme_minimal(base_size = 12)
save_plot(p1, "01_barcode_sharing_histogram", w = 7, h = 4.5)

# ---- 2. BARCODES PER SAMPLE, FACETED BY MOUSE --------------------------
per_sample <- long |>
  dplyr::distinct(mouse, sample_type, sample, barcode) |>
  dplyr::count(mouse, sample_type, sample, name = "n_barcodes") |>
  dplyr::mutate(sample_type = factor(sample_type, levels = sample_type_levels))

readr::write_csv(per_sample, file.path(output_dir, "barcodes_per_sample.csv"))

p2 <- per_sample |>
  ggplot(aes(x = n_barcodes, y = sample_type, fill = factor(mouse))) +
  geom_col() +
  geom_text(aes(label = n_barcodes), hjust = -0.15, size = 3) +
  facet_wrap(~ paste("Mouse", mouse), nrow = 1) +
  scale_y_discrete(drop = FALSE, limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "# barcodes", y = NULL,
       title = "Barcodes per sample, per mouse") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
save_plot(p2, "02_barcodes_per_sample_per_mouse",
          w = max(10, 3.2 * dplyr::n_distinct(meta$mouse)), h = 5)

# ---- 3a. UPSET PER MOUSE (ALL SAMPLES) ---------------------------------
for (m in sort(unique(meta$mouse))) {
  m_meta    <- meta |> dplyr::filter(mouse == m)
  m_samples <- m_meta$sample
  rename_to <- m_meta$sample_type
  if (anyDuplicated(rename_to)) rename_to <- m_samples

  pa <- wide |>
    dplyr::select(barcode, all_of(m_samples)) |>
    dplyr::rename_with(~ rename_to, all_of(m_samples)) |>
    dplyr::mutate(dplyr::across(all_of(rename_to), ~ .x > 0)) |>
    dplyr::filter(rowSums(dplyr::across(all_of(rename_to))) > 0)

  if (length(rename_to) < 2) {
    message("Mouse ", m, ": only one sample, skipping UpSet.")
    next
  }

  p_upset <- ComplexUpset::upset(
    pa, intersect = rename_to,
    name = "Sample",
    width_ratio = 0.25,
    sort_intersections_by = "cardinality",
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(text = list(size = 3))
    )
  ) +
    patchwork::plot_annotation(
      title = paste0("Mouse ", m, " - barcode intersections across all samples")
    )

  save_plot(p_upset, sprintf("03a_upset_mouse%d", m),
            w = max(8, 1.1 * length(rename_to) + 5), h = 5.5)
}

# ---- 3b. VENN PER MOUSE x POPULATION -----------------------------------
# One Venn per (mouse, population). Sets within = bone/source samples.
for (m in sort(unique(meta$mouse))) {
  for (pop in populations) {
    sel <- meta |>
      dplyr::filter(mouse == m, population == pop) |>
      dplyr::arrange(source)

    if (nrow(sel) < 2) {
      message("Mouse ", m, " / ", pop, ": ", nrow(sel),
              " sample(s) - skipping Venn.")
      next
    }
    if (nrow(sel) > 7) {
      message("Mouse ", m, " / ", pop, ": ", nrow(sel),
              " samples - too many for a clean Venn (>7); skipping. ",
              "See UpSet plot 03a for this mouse.")
      next
    }

    set_labels <- sel$source
    if (anyDuplicated(set_labels)) set_labels <- sel$sample

    bc_sets <- setNames(
      lapply(sel$sample, function(s) wide$barcode[wide[[s]] > 0]),
      set_labels
    )

    p_venn <- ggVennDiagram::ggVennDiagram(
        bc_sets, label_alpha = 0, set_size = 3.5, label = "count"
      ) +
      scale_fill_distiller(palette = "Blues", direction = 1) +
      scale_x_continuous(expand = expansion(mult = 0.2)) +
      labs(title = sprintf("Mouse %d - %s - barcode overlap across bones",
                           m, pop)) +
      theme(legend.position = "none")

    save_plot(p_venn, sprintf("03b_venn_mouse%d_%s", m, pop),
              w = 6.5, h = 5.5)
  }
}

# ---- 4. HEATMAP: barcodes per mouse x sample_type ----------------------
heat <- per_sample |>
  dplyr::group_by(mouse, sample_type) |>
  dplyr::summarise(n_barcodes = sum(n_barcodes), .groups = "drop")

heat_wide <- heat |>
  tidyr::pivot_wider(names_from = sample_type, values_from = n_barcodes, values_fill = 0)
readr::write_csv(heat_wide, file.path(output_dir, "heatmap_barcodes_mouse_x_sampletype.csv"))

text_thresh <- 0.5 * max(heat$n_barcodes, na.rm = TRUE)
p4 <- heat |>
  ggplot(aes(x = sample_type, y = factor(mouse), fill = n_barcodes)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_barcodes,
                color = n_barcodes > text_thresh),
            size = 3.5) +
  scale_color_manual(values = c("FALSE" = "white", "TRUE" = "black"), guide = "none") +
  scale_fill_viridis_c(option = "viridis", name = "# barcodes") +
  labs(x = NULL, y = "Mouse",
       title = "Number of barcodes per sample type, all mice") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
save_plot(p4, "04_heatmap_all_mice",
          w = max(7, 1.1 * length(sample_type_levels) + 3), h = 4.5)

message("All outputs written to: ", output_dir)
