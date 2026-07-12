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
input_dir   <- "C:/Users/fc809/Downloads/T_region_result/results_T_regions"                              # change to your folder
output_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"

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


#####SUMMARY_LOOK_VENN_HEATMAP_########


# barcode_visualizations.R
# Reads combined_barcodes_wide.tsv and produces:
#   1.  Histogram of barcode sharing across samples
#   2.  Barcodes per sample, faceted by mouse
#   3a. UpSet plot per mouse (all samples)
#   3b. Venn diagram per mouse PER hematopoietic POPULATION
#       (sets within each Venn = the different bone/source samples)
#   4.  Heatmap of barcode counts (mouse x sample_type)
# Each plot saved as PNG and PDF; summary tables saved as CSV.
#install.packages("ComplexUpset")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)   # install.packages("ggVennDiagram")
  library(ComplexUpset)    # install.packages("ComplexUpset")
  library(scales)
})

# ---- USER SETTINGS ------------------------------------------------------
input_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"   # in working dir, or give full path
output_dir <- "C:/Users/fc809/Downloads/T_region_result/"

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

#########distance_analysis#################3
# thymus_distance.R
# Per mouse, compute Euclidean and Bray-Curtis distance from each
# (population x bone) source sample to that mouse's Thymus_DN1 sample,
# using barcode frequencies (UMI counts normalized to proportions).
# Aggregate across mice, run Friedman test, and produce:
#   05a  ranked dot+bar plot of all sources (one per metric)
#   05b  faceted-by-population plot (one per metric)
# plus CSVs of the per-mouse distances and the Friedman results.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- USER SETTINGS ------------------------------------------------------
input_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"
output_dir <- "C:/Users/fc809/Downloads/T_region_result/"

# Destination sample: matched against sample_type by str_detect (so e.g.
# 'Thymus_DN1' will catch both 'Thymus_DN1' and 'Thymus_DN1_CA').
destination_pattern <- "Thymus_DN1"

# Source samples kept: must have one of these populations AND one of these bones.
populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")
# -------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# ---- LOAD & PARSE -------------------------------------------------------
wide <- readr::read_tsv(input_file, show_col_types = FALSE)
summary_cols <- c("barcode", "n_samples_detected", "total_umi")
sample_cols  <- setdiff(colnames(wide), summary_cols)

# Sample naming: {mouse}_{population}_{bone}_{condition}_S{##}
# We parse by splitting on underscores and taking specific tokens, which avoids
# the previous bug where "_CA" condition leaked into the bone field.
meta <- tibble::tibble(sample = sample_cols) |>
  dplyr::mutate(
    parts       = stringr::str_split(sample, "_"),
    mouse       = suppressWarnings(as.integer(purrr::map_chr(parts, 1L))),
    population  = purrr::map_chr(parts, ~ if (length(.x) >= 2) .x[2] else NA_character_),
    bone        = purrr::map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
    condition   = purrr::map_chr(parts, ~ if (length(.x) >= 4) .x[4] else NA_character_),
    sample_type = sample |>
      stringr::str_remove("^[^_]+_") |>   # strip whatever-came-first (digits or "WT")
      stringr::str_remove("_S\\d+$")       # strip trailing _S##
  ) |>
  dplyr::select(-parts)

# Destination: filter to non-NA mouse so WT controls (mouse = NA) are excluded
dest_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                stringr::str_detect(sample_type, destination_pattern))

src_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                population %in% populations,
                bone        %in% bones)

message("Destination samples (", destination_pattern, "):")
print(dest_meta |> dplyr::select(mouse, sample))
message("Source samples kept: ", nrow(src_meta),
        " across mice ", paste(sort(unique(src_meta$mouse)), collapse = ", "))
message("Populations in sources: ",
        paste(sort(unique(src_meta$population)), collapse = ", "))
message("Bones in sources: ",
        paste(sort(unique(src_meta$bone)), collapse = ", "))

if (anyDuplicated(dest_meta$mouse)) {
  warning("Multiple destination samples found in some mice; the first per mouse will be used.")
  dest_meta <- dest_meta |> dplyr::group_by(mouse) |> dplyr::slice(1) |> dplyr::ungroup()
}

# ---- DISTANCE FUNCTIONS -------------------------------------------------
euclidean   <- function(a, b) sqrt(sum((a - b)^2))
bray_curtis <- function(a, b) {
  denom <- sum(a + b)
  if (denom == 0) NA_real_ else sum(abs(a - b)) / denom
}

# ---- PER-MOUSE COMPUTATION ----------------------------------------------
results <- list()
for (m in sort(unique(dest_meta$mouse))) {
  dest_col <- dest_meta$sample[dest_meta$mouse == m]
  m_src    <- src_meta |> dplyr::filter(mouse == m)
  if (nrow(m_src) == 0) {
    message("Mouse ", m, ": destination found but no qualifying source samples; skipping.")
    next
  }
  m_cols <- c(dest_col, m_src$sample)
  
  # Submatrix for this mouse, restricted to barcodes detected in any of these samples
  M <- as.matrix(wide[, m_cols])
  M <- M[rowSums(M) > 0, , drop = FALSE]
  
  # Normalize each column to sum = 1 (frequencies)
  totals <- colSums(M); totals[totals == 0] <- 1
  Fmat <- sweep(M, 2, totals, "/")
  
  dest_freq <- Fmat[, dest_col]
  for (i in seq_len(nrow(m_src))) {
    src_col  <- m_src$sample[i]
    src_freq <- Fmat[, src_col]
    results[[length(results) + 1]] <- tibble::tibble(
      mouse       = m,
      population  = m_src$population[i],
      bone        = m_src$bone[i],
      source      = paste(m_src$population[i], m_src$bone[i], sep = "_"),
      sample      = src_col,
      euclidean   = euclidean(src_freq, dest_freq),
      bray_curtis = bray_curtis(src_freq, dest_freq)
    )
  }
}

dist_df <- dplyr::bind_rows(results)
if (nrow(dist_df) == 0) stop("No distances computed; check destination_pattern and sample naming.")

# Collapse any duplicate (mouse, source) entries by averaging
dist_df <- dist_df |>
  dplyr::group_by(mouse, population, bone, source) |>
  dplyr::summarise(euclidean   = mean(euclidean),
                   bray_curtis = mean(bray_curtis),
                   .groups     = "drop")

readr::write_csv(dist_df, file.path(output_dir, "thymus_distances_per_mouse.csv"))

# ---- FRIEDMAN TEST (complete blocks only) -------------------------------
friedman_one <- function(metric) {
  d <- dist_df |> dplyr::select(mouse, source, value = all_of(metric))
  n_mice <- dplyr::n_distinct(d$mouse)
  complete <- d |> dplyr::count(source) |>
    dplyr::filter(n == n_mice) |> dplyr::pull(source)
  if (length(complete) < 2) {
    return(tibble::tibble(metric = metric, n_sources_tested = length(complete),
                          n_mice = n_mice, statistic = NA_real_,
                          df = NA_real_, p_value = NA_real_))
  }
  d2 <- d |> dplyr::filter(source %in% complete) |>
    dplyr::mutate(mouse = factor(mouse), source = factor(source))
  ft <- stats::friedman.test(value ~ source | mouse, data = d2)
  tibble::tibble(metric = metric, n_sources_tested = length(complete),
                 n_mice = n_mice,
                 statistic = unname(ft$statistic),
                 df = unname(ft$parameter),
                 p_value = ft$p.value)
}
friedman_df <- dplyr::bind_rows(lapply(c("euclidean", "bray_curtis"), friedman_one))
readr::write_csv(friedman_df, file.path(output_dir, "thymus_distance_friedman.csv"))
print(friedman_df)

# ---- PLOTS --------------------------------------------------------------
pop_palette <- c(LSK = "#2A9D8F", MPP4 = "#E9C46A",
                 CLP = "#8D6CAB", GMP = "#F4A261")

metric_label <- function(mtr) {
  if (mtr == "euclidean")
    "Euclidean distance to Thymus_DN1\n(lower = more similar)"
  else
    "Bray-Curtis dissimilarity to Thymus_DN1\n(lower = more similar)"
}

# --- Plot A: ranked, single panel ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |> dplyr::rename(value = all_of(mtr))
  src_order <- d |> dplyr::group_by(source) |>
    dplyr::summarise(m = mean(value), .groups = "drop") |>
    dplyr::arrange(m) |> dplyr::pull(source)
  
  agg <- d |> dplyr::group_by(source, population) |>
    dplyr::summarise(mean_d = mean(value),
                     sem    = sd(value) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(source = factor(source, levels = src_order))
  
  ft <- friedman_df |> dplyr::filter(metric == mtr)
  subtitle <- if (!is.na(ft$p_value))
    sprintf("Friedman across %d complete-block sources (n=%d mice): chi2 = %.2f, p = %.3g",
            ft$n_sources_tested, ft$n_mice, ft$statistic, ft$p_value)
  else "Friedman not run (fewer than 2 complete-block sources)"
  
  d_pts <- d |> dplyr::mutate(source = factor(source, levels = src_order))
  
  p <- ggplot(agg, aes(x = source, y = mean_d, fill = population)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.78) +
    geom_errorbar(aes(ymin = mean_d - sem, ymax = mean_d + sem),
                  width = 0.3, linewidth = 0.4) +
    geom_point(data = d_pts, inherit.aes = FALSE,
               mapping = aes(x = source, y = value, color = factor(mouse)),
               size = 1.9, alpha = 0.9,
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    scale_fill_manual(values = pop_palette, name = "Population") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Which (population x bone) source is closest to Thymus_DN1?",
         subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05a_thymus_distance_ranked_", mtr),
            w = max(11, 0.55 * length(src_order) + 3), h = 5.8)
}

# --- Plot B: faceted by population, lines connecting mice across bones ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |>
    dplyr::rename(value = all_of(mtr)) |>
    dplyr::mutate(bone = factor(bone, levels = bones))
  
  p <- d |>
    ggplot(aes(x = bone, y = value)) +
    stat_summary(aes(fill = population), fun = mean, geom = "col",
                 color = "black", linewidth = 0.3, alpha = 0.75) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.3, linewidth = 0.4) +
    geom_line(aes(group = mouse, color = factor(mouse)), alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(color = factor(mouse)), size = 2, alpha = 0.9) +
    facet_wrap(~ population, nrow = 1) +
    scale_fill_manual(values = pop_palette, guide = "none") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Distance to Thymus_DN1 by bone, faceted by population") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05b_thymus_distance_faceted_", mtr),
            w = max(13, 3.4 * length(populations)), h = 5.2)
}

message("Done. Outputs in: ", output_dir)

###############################Blood_T_Cells################################

# ---- USER SETTINGS ------------------------------------------------------
input_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"
output_dir <- "C:/Users/fc809/Downloads/T_region_result/"

# Destination sample: matched against sample_type by str_detect (so e.g.
# 'Blood_T' will catch both 'Blood_T').
destination_pattern <- "Blood_T"

# Source samples kept: must have one of these populations AND one of these bones.
populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")
# -------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# ---- LOAD & PARSE -------------------------------------------------------
wide <- readr::read_tsv(input_file, show_col_types = FALSE)
summary_cols <- c("barcode", "n_samples_detected", "total_umi")
sample_cols  <- setdiff(colnames(wide), summary_cols)

# Sample naming: {mouse}_{population}_{bone}_{condition}_S{##}
# We parse by splitting on underscores and taking specific tokens, which avoids
# the previous bug where "_CA" condition leaked into the bone field.
meta <- tibble::tibble(sample = sample_cols) |>
  dplyr::mutate(
    parts       = stringr::str_split(sample, "_"),
    mouse       = suppressWarnings(as.integer(purrr::map_chr(parts, 1L))),
    population  = purrr::map_chr(parts, ~ if (length(.x) >= 2) .x[2] else NA_character_),
    bone        = purrr::map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
    condition   = purrr::map_chr(parts, ~ if (length(.x) >= 4) .x[4] else NA_character_),
    sample_type = sample |>
      stringr::str_remove("^[^_]+_") |>   # strip whatever-came-first (digits or "WT")
      stringr::str_remove("_S\\d+$")       # strip trailing _S##
  ) |>
  dplyr::select(-parts)

# Destination: filter to non-NA mouse so WT controls (mouse = NA) are excluded
dest_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                stringr::str_detect(sample_type, destination_pattern))

src_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                population %in% populations,
                bone        %in% bones)

message("Destination samples (", destination_pattern, "):")
print(dest_meta |> dplyr::select(mouse, sample))
message("Source samples kept: ", nrow(src_meta),
        " across mice ", paste(sort(unique(src_meta$mouse)), collapse = ", "))
message("Populations in sources: ",
        paste(sort(unique(src_meta$population)), collapse = ", "))
message("Bones in sources: ",
        paste(sort(unique(src_meta$bone)), collapse = ", "))

if (anyDuplicated(dest_meta$mouse)) {
  warning("Multiple destination samples found in some mice; the first per mouse will be used.")
  dest_meta <- dest_meta |> dplyr::group_by(mouse) |> dplyr::slice(1) |> dplyr::ungroup()
}

# ---- DISTANCE FUNCTIONS -------------------------------------------------
euclidean   <- function(a, b) sqrt(sum((a - b)^2))
bray_curtis <- function(a, b) {
  denom <- sum(a + b)
  if (denom == 0) NA_real_ else sum(abs(a - b)) / denom
}

# ---- PER-MOUSE COMPUTATION ----------------------------------------------
results <- list()
for (m in sort(unique(dest_meta$mouse))) {
  dest_col <- dest_meta$sample[dest_meta$mouse == m]
  m_src    <- src_meta |> dplyr::filter(mouse == m)
  if (nrow(m_src) == 0) {
    message("Mouse ", m, ": destination found but no qualifying source samples; skipping.")
    next
  }
  m_cols <- c(dest_col, m_src$sample)
  
  # Submatrix for this mouse, restricted to barcodes detected in any of these samples
  M <- as.matrix(wide[, m_cols])
  M <- M[rowSums(M) > 0, , drop = FALSE]
  
  # Normalize each column to sum = 1 (frequencies)
  totals <- colSums(M); totals[totals == 0] <- 1
  Fmat <- sweep(M, 2, totals, "/")
  
  dest_freq <- Fmat[, dest_col]
  for (i in seq_len(nrow(m_src))) {
    src_col  <- m_src$sample[i]
    src_freq <- Fmat[, src_col]
    results[[length(results) + 1]] <- tibble::tibble(
      mouse       = m,
      population  = m_src$population[i],
      bone        = m_src$bone[i],
      source      = paste(m_src$population[i], m_src$bone[i], sep = "_"),
      sample      = src_col,
      euclidean   = euclidean(src_freq, dest_freq),
      bray_curtis = bray_curtis(src_freq, dest_freq)
    )
  }
}

dist_df <- dplyr::bind_rows(results)
if (nrow(dist_df) == 0) stop("No distances computed; check destination_pattern and sample naming.")

# Collapse any duplicate (mouse, source) entries by averaging
dist_df <- dist_df |>
  dplyr::group_by(mouse, population, bone, source) |>
  dplyr::summarise(euclidean   = mean(euclidean),
                   bray_curtis = mean(bray_curtis),
                   .groups     = "drop")

readr::write_csv(dist_df, file.path(output_dir, "blood_T_distances_per_mouse.csv"))

# ---- FRIEDMAN TEST (complete blocks only) -------------------------------
friedman_one <- function(metric) {
  d <- dist_df |> dplyr::select(mouse, source, value = all_of(metric))
  n_mice <- dplyr::n_distinct(d$mouse)
  complete <- d |> dplyr::count(source) |>
    dplyr::filter(n == n_mice) |> dplyr::pull(source)
  if (length(complete) < 2) {
    return(tibble::tibble(metric = metric, n_sources_tested = length(complete),
                          n_mice = n_mice, statistic = NA_real_,
                          df = NA_real_, p_value = NA_real_))
  }
  d2 <- d |> dplyr::filter(source %in% complete) |>
    dplyr::mutate(mouse = factor(mouse), source = factor(source))
  ft <- stats::friedman.test(value ~ source | mouse, data = d2)
  tibble::tibble(metric = metric, n_sources_tested = length(complete),
                 n_mice = n_mice,
                 statistic = unname(ft$statistic),
                 df = unname(ft$parameter),
                 p_value = ft$p.value)
}
friedman_df <- dplyr::bind_rows(lapply(c("euclidean", "bray_curtis"), friedman_one))
readr::write_csv(friedman_df, file.path(output_dir, "blood_T_distance_friedman.csv"))
print(friedman_df)

# ---- PLOTS --------------------------------------------------------------
pop_palette <- c(LSK = "#2A9D8F", MPP4 = "#E9C46A",
                 CLP = "#8D6CAB", GMP = "#F4A261")

metric_label <- function(mtr) {
  if (mtr == "euclidean")
    "Euclidean distance to Blood_T\n(lower = more similar)"
  else
    "Bray-Curtis dissimilarity to Blood_T\n(lower = more similar)"
}

# --- Plot A: ranked, single panel ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |> dplyr::rename(value = all_of(mtr))
  src_order <- d |> dplyr::group_by(source) |>
    dplyr::summarise(m = mean(value), .groups = "drop") |>
    dplyr::arrange(m) |> dplyr::pull(source)
  
  agg <- d |> dplyr::group_by(source, population) |>
    dplyr::summarise(mean_d = mean(value),
                     sem    = sd(value) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(source = factor(source, levels = src_order))
  
  ft <- friedman_df |> dplyr::filter(metric == mtr)
  subtitle <- if (!is.na(ft$p_value))
    sprintf("Friedman across %d complete-block sources (n=%d mice): chi2 = %.2f, p = %.3g",
            ft$n_sources_tested, ft$n_mice, ft$statistic, ft$p_value)
  else "Friedman not run (fewer than 2 complete-block sources)"
  
  d_pts <- d |> dplyr::mutate(source = factor(source, levels = src_order))
  
  p <- ggplot(agg, aes(x = source, y = mean_d, fill = population)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.78) +
    geom_errorbar(aes(ymin = mean_d - sem, ymax = mean_d + sem),
                  width = 0.3, linewidth = 0.4) +
    geom_point(data = d_pts, inherit.aes = FALSE,
               mapping = aes(x = source, y = value, color = factor(mouse)),
               size = 1.9, alpha = 0.9,
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    scale_fill_manual(values = pop_palette, name = "Population") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Which (population x bone) source is closest to Blood_T?",
         subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05a_Blood_T_distance_ranked_", mtr),
            w = max(11, 0.55 * length(src_order) + 3), h = 5.8)
}

# --- Plot B: faceted by population, lines connecting mice across bones ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |>
    dplyr::rename(value = all_of(mtr)) |>
    dplyr::mutate(bone = factor(bone, levels = bones))
  
  p <- d |>
    ggplot(aes(x = bone, y = value)) +
    stat_summary(aes(fill = population), fun = mean, geom = "col",
                 color = "black", linewidth = 0.3, alpha = 0.75) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.3, linewidth = 0.4) +
    geom_line(aes(group = mouse, color = factor(mouse)), alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(color = factor(mouse)), size = 2, alpha = 0.9) +
    facet_wrap(~ population, nrow = 1) +
    scale_fill_manual(values = pop_palette, guide = "none") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Distance to Blood_T by bone, faceted by population") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05b_Blood_T_distance_faceted_", mtr),
            w = max(13, 3.4 * length(populations)), h = 5.2)
}

message("Done. Outputs in: ", output_dir)

###############################Blood_Granulocytes_Cells################################

# ---- USER SETTINGS ------------------------------------------------------
input_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"
output_dir <- "C:/Users/fc809/Downloads/T_region_result/"

# Destination sample: matched against sample_type by str_detect (so e.g.
# 'Blood_Granulocytes' will catch both 'Blood_Granulocytes').
destination_pattern <- "Blood_Granulocytes"

# Source samples kept: must have one of these populations AND one of these bones.
populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")
# -------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# ---- LOAD & PARSE -------------------------------------------------------
wide <- readr::read_tsv(input_file, show_col_types = FALSE)
summary_cols <- c("barcode", "n_samples_detected", "total_umi")
sample_cols  <- setdiff(colnames(wide), summary_cols)

# Sample naming: {mouse}_{population}_{bone}_{condition}_S{##}
# We parse by splitting on underscores and taking specific tokens, which avoids
# the previous bug where "_CA" condition leaked into the bone field.
meta <- tibble::tibble(sample = sample_cols) |>
  dplyr::mutate(
    parts       = stringr::str_split(sample, "_"),
    mouse       = suppressWarnings(as.integer(purrr::map_chr(parts, 1L))),
    population  = purrr::map_chr(parts, ~ if (length(.x) >= 2) .x[2] else NA_character_),
    bone        = purrr::map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
    condition   = purrr::map_chr(parts, ~ if (length(.x) >= 4) .x[4] else NA_character_),
    sample_type = sample |>
      stringr::str_remove("^[^_]+_") |>   # strip whatever-came-first (digits or "WT")
      stringr::str_remove("_S\\d+$")       # strip trailing _S##
  ) |>
  dplyr::select(-parts)

# Destination: filter to non-NA mouse so WT controls (mouse = NA) are excluded
dest_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                stringr::str_detect(sample_type, destination_pattern))

src_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                population %in% populations,
                bone        %in% bones)

message("Destination samples (", destination_pattern, "):")
print(dest_meta |> dplyr::select(mouse, sample))
message("Source samples kept: ", nrow(src_meta),
        " across mice ", paste(sort(unique(src_meta$mouse)), collapse = ", "))
message("Populations in sources: ",
        paste(sort(unique(src_meta$population)), collapse = ", "))
message("Bones in sources: ",
        paste(sort(unique(src_meta$bone)), collapse = ", "))

if (anyDuplicated(dest_meta$mouse)) {
  warning("Multiple destination samples found in some mice; the first per mouse will be used.")
  dest_meta <- dest_meta |> dplyr::group_by(mouse) |> dplyr::slice(1) |> dplyr::ungroup()
}

# ---- DISTANCE FUNCTIONS -------------------------------------------------
euclidean   <- function(a, b) sqrt(sum((a - b)^2))
bray_curtis <- function(a, b) {
  denom <- sum(a + b)
  if (denom == 0) NA_real_ else sum(abs(a - b)) / denom
}

# ---- PER-MOUSE COMPUTATION ----------------------------------------------
results <- list()
for (m in sort(unique(dest_meta$mouse))) {
  dest_col <- dest_meta$sample[dest_meta$mouse == m]
  m_src    <- src_meta |> dplyr::filter(mouse == m)
  if (nrow(m_src) == 0) {
    message("Mouse ", m, ": destination found but no qualifying source samples; skipping.")
    next
  }
  m_cols <- c(dest_col, m_src$sample)
  
  # Submatrix for this mouse, restricted to barcodes detected in any of these samples
  M <- as.matrix(wide[, m_cols])
  M <- M[rowSums(M) > 0, , drop = FALSE]
  
  # Normalize each column to sum = 1 (frequencies)
  totals <- colSums(M); totals[totals == 0] <- 1
  Fmat <- sweep(M, 2, totals, "/")
  
  dest_freq <- Fmat[, dest_col]
  for (i in seq_len(nrow(m_src))) {
    src_col  <- m_src$sample[i]
    src_freq <- Fmat[, src_col]
    results[[length(results) + 1]] <- tibble::tibble(
      mouse       = m,
      population  = m_src$population[i],
      bone        = m_src$bone[i],
      source      = paste(m_src$population[i], m_src$bone[i], sep = "_"),
      sample      = src_col,
      euclidean   = euclidean(src_freq, dest_freq),
      bray_curtis = bray_curtis(src_freq, dest_freq)
    )
  }
}

dist_df <- dplyr::bind_rows(results)
if (nrow(dist_df) == 0) stop("No distances computed; check destination_pattern and sample naming.")

# Collapse any duplicate (mouse, source) entries by averaging
dist_df <- dist_df |>
  dplyr::group_by(mouse, population, bone, source) |>
  dplyr::summarise(euclidean   = mean(euclidean),
                   bray_curtis = mean(bray_curtis),
                   .groups     = "drop")

readr::write_csv(dist_df, file.path(output_dir, "Blood_Granulocytes_distances_per_mouse.csv"))

# ---- FRIEDMAN TEST (complete blocks only) -------------------------------
friedman_one <- function(metric) {
  d <- dist_df |> dplyr::select(mouse, source, value = all_of(metric))
  n_mice <- dplyr::n_distinct(d$mouse)
  complete <- d |> dplyr::count(source) |>
    dplyr::filter(n == n_mice) |> dplyr::pull(source)
  if (length(complete) < 2) {
    return(tibble::tibble(metric = metric, n_sources_tested = length(complete),
                          n_mice = n_mice, statistic = NA_real_,
                          df = NA_real_, p_value = NA_real_))
  }
  d2 <- d |> dplyr::filter(source %in% complete) |>
    dplyr::mutate(mouse = factor(mouse), source = factor(source))
  ft <- stats::friedman.test(value ~ source | mouse, data = d2)
  tibble::tibble(metric = metric, n_sources_tested = length(complete),
                 n_mice = n_mice,
                 statistic = unname(ft$statistic),
                 df = unname(ft$parameter),
                 p_value = ft$p.value)
}
friedman_df <- dplyr::bind_rows(lapply(c("euclidean", "bray_curtis"), friedman_one))
readr::write_csv(friedman_df, file.path(output_dir, "Blood_Granulocytes_distance_friedman.csv"))
print(friedman_df)

# ---- PLOTS --------------------------------------------------------------
pop_palette <- c(LSK = "#2A9D8F", MPP4 = "#E9C46A",
                 CLP = "#8D6CAB", GMP = "#F4A261")

metric_label <- function(mtr) {
  if (mtr == "euclidean")
    "Euclidean distance to Blood_Granulocytes\n(lower = more similar)"
  else
    "Bray-Curtis dissimilarity to Blood_Granulocytes\n(lower = more similar)"
}

# --- Plot A: ranked, single panel ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |> dplyr::rename(value = all_of(mtr))
  src_order <- d |> dplyr::group_by(source) |>
    dplyr::summarise(m = mean(value), .groups = "drop") |>
    dplyr::arrange(m) |> dplyr::pull(source)
  
  agg <- d |> dplyr::group_by(source, population) |>
    dplyr::summarise(mean_d = mean(value),
                     sem    = sd(value) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(source = factor(source, levels = src_order))
  
  ft <- friedman_df |> dplyr::filter(metric == mtr)
  subtitle <- if (!is.na(ft$p_value))
    sprintf("Friedman across %d complete-block sources (n=%d mice): chi2 = %.2f, p = %.3g",
            ft$n_sources_tested, ft$n_mice, ft$statistic, ft$p_value)
  else "Friedman not run (fewer than 2 complete-block sources)"
  
  d_pts <- d |> dplyr::mutate(source = factor(source, levels = src_order))
  
  p <- ggplot(agg, aes(x = source, y = mean_d, fill = population)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.78) +
    geom_errorbar(aes(ymin = mean_d - sem, ymax = mean_d + sem),
                  width = 0.3, linewidth = 0.4) +
    geom_point(data = d_pts, inherit.aes = FALSE,
               mapping = aes(x = source, y = value, color = factor(mouse)),
               size = 1.9, alpha = 0.9,
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    scale_fill_manual(values = pop_palette, name = "Population") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Which (population x bone) source is closest to Blood_Granulocytes?",
         subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05a_Blood_Granulocytes_distance_ranked_", mtr),
            w = max(11, 0.55 * length(src_order) + 3), h = 5.8)
}

# --- Plot B: faceted by population, lines connecting mice across bones ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |>
    dplyr::rename(value = all_of(mtr)) |>
    dplyr::mutate(bone = factor(bone, levels = bones))
  
  p <- d |>
    ggplot(aes(x = bone, y = value)) +
    stat_summary(aes(fill = population), fun = mean, geom = "col",
                 color = "black", linewidth = 0.3, alpha = 0.75) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.3, linewidth = 0.4) +
    geom_line(aes(group = mouse, color = factor(mouse)), alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(color = factor(mouse)), size = 2, alpha = 0.9) +
    facet_wrap(~ population, nrow = 1) +
    scale_fill_manual(values = pop_palette, guide = "none") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Distance to Blood_Granulocytes by bone, faceted by population") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05b_Blood_Granulocytes_distance_faceted_", mtr),
            w = max(13, 3.4 * length(populations)), h = 5.2)
}

message("Done. Outputs in: ", output_dir)

###############################Heart_Granulocytes_Cells################################

# ---- USER SETTINGS ------------------------------------------------------
input_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"
output_dir <- "C:/Users/fc809/Downloads/T_region_result/"

# Destination sample: matched against sample_type by str_detect (so e.g.
# 'Heart_Granulocytes' will catch both 'Heart_Granulocytes').
destination_pattern <- "Heart_Granulocytes"

# Source samples kept: must have one of these populations AND one of these bones.
populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")
# -------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# ---- LOAD & PARSE -------------------------------------------------------
wide <- readr::read_tsv(input_file, show_col_types = FALSE)
summary_cols <- c("barcode", "n_samples_detected", "total_umi")
sample_cols  <- setdiff(colnames(wide), summary_cols)

# Sample naming: {mouse}_{population}_{bone}_{condition}_S{##}
# We parse by splitting on underscores and taking specific tokens, which avoids
# the previous bug where "_CA" condition leaked into the bone field.
meta <- tibble::tibble(sample = sample_cols) |>
  dplyr::mutate(
    parts       = stringr::str_split(sample, "_"),
    mouse       = suppressWarnings(as.integer(purrr::map_chr(parts, 1L))),
    population  = purrr::map_chr(parts, ~ if (length(.x) >= 2) .x[2] else NA_character_),
    bone        = purrr::map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
    condition   = purrr::map_chr(parts, ~ if (length(.x) >= 4) .x[4] else NA_character_),
    sample_type = sample |>
      stringr::str_remove("^[^_]+_") |>   # strip whatever-came-first (digits or "WT")
      stringr::str_remove("_S\\d+$")       # strip trailing _S##
  ) |>
  dplyr::select(-parts)

# Destination: filter to non-NA mouse so WT controls (mouse = NA) are excluded
dest_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                stringr::str_detect(sample_type, destination_pattern))

src_meta <- meta |>
  dplyr::filter(!is.na(mouse),
                population %in% populations,
                bone        %in% bones)

message("Destination samples (", destination_pattern, "):")
print(dest_meta |> dplyr::select(mouse, sample))
message("Source samples kept: ", nrow(src_meta),
        " across mice ", paste(sort(unique(src_meta$mouse)), collapse = ", "))
message("Populations in sources: ",
        paste(sort(unique(src_meta$population)), collapse = ", "))
message("Bones in sources: ",
        paste(sort(unique(src_meta$bone)), collapse = ", "))

if (anyDuplicated(dest_meta$mouse)) {
  warning("Multiple destination samples found in some mice; the first per mouse will be used.")
  dest_meta <- dest_meta |> dplyr::group_by(mouse) |> dplyr::slice(1) |> dplyr::ungroup()
}

# ---- DISTANCE FUNCTIONS -------------------------------------------------
euclidean   <- function(a, b) sqrt(sum((a - b)^2))
bray_curtis <- function(a, b) {
  denom <- sum(a + b)
  if (denom == 0) NA_real_ else sum(abs(a - b)) / denom
}

# ---- PER-MOUSE COMPUTATION ----------------------------------------------
results <- list()
for (m in sort(unique(dest_meta$mouse))) {
  dest_col <- dest_meta$sample[dest_meta$mouse == m]
  m_src    <- src_meta |> dplyr::filter(mouse == m)
  if (nrow(m_src) == 0) {
    message("Mouse ", m, ": destination found but no qualifying source samples; skipping.")
    next
  }
  m_cols <- c(dest_col, m_src$sample)
  
  # Submatrix for this mouse, restricted to barcodes detected in any of these samples
  M <- as.matrix(wide[, m_cols])
  M <- M[rowSums(M) > 0, , drop = FALSE]
  
  # Normalize each column to sum = 1 (frequencies)
  totals <- colSums(M); totals[totals == 0] <- 1
  Fmat <- sweep(M, 2, totals, "/")
  
  dest_freq <- Fmat[, dest_col]
  for (i in seq_len(nrow(m_src))) {
    src_col  <- m_src$sample[i]
    src_freq <- Fmat[, src_col]
    results[[length(results) + 1]] <- tibble::tibble(
      mouse       = m,
      population  = m_src$population[i],
      bone        = m_src$bone[i],
      source      = paste(m_src$population[i], m_src$bone[i], sep = "_"),
      sample      = src_col,
      euclidean   = euclidean(src_freq, dest_freq),
      bray_curtis = bray_curtis(src_freq, dest_freq)
    )
  }
}

dist_df <- dplyr::bind_rows(results)
if (nrow(dist_df) == 0) stop("No distances computed; check destination_pattern and sample naming.")

# Collapse any duplicate (mouse, source) entries by averaging
dist_df <- dist_df |>
  dplyr::group_by(mouse, population, bone, source) |>
  dplyr::summarise(euclidean   = mean(euclidean),
                   bray_curtis = mean(bray_curtis),
                   .groups     = "drop")

readr::write_csv(dist_df, file.path(output_dir, "Heart_Granulocytes_distances_per_mouse.csv"))

# ---- FRIEDMAN TEST (complete blocks only) -------------------------------
friedman_one <- function(metric) {
  d <- dist_df |> dplyr::select(mouse, source, value = all_of(metric))
  n_mice <- dplyr::n_distinct(d$mouse)
  complete <- d |> dplyr::count(source) |>
    dplyr::filter(n == n_mice) |> dplyr::pull(source)
  if (length(complete) < 2) {
    return(tibble::tibble(metric = metric, n_sources_tested = length(complete),
                          n_mice = n_mice, statistic = NA_real_,
                          df = NA_real_, p_value = NA_real_))
  }
  d2 <- d |> dplyr::filter(source %in% complete) |>
    dplyr::mutate(mouse = factor(mouse), source = factor(source))
  ft <- stats::friedman.test(value ~ source | mouse, data = d2)
  tibble::tibble(metric = metric, n_sources_tested = length(complete),
                 n_mice = n_mice,
                 statistic = unname(ft$statistic),
                 df = unname(ft$parameter),
                 p_value = ft$p.value)
}
friedman_df <- dplyr::bind_rows(lapply(c("euclidean", "bray_curtis"), friedman_one))
readr::write_csv(friedman_df, file.path(output_dir, "Heart_Granulocytes_distance_friedman.csv"))
print(friedman_df)

# ---- PLOTS --------------------------------------------------------------
pop_palette <- c(LSK = "#2A9D8F", MPP4 = "#E9C46A",
                 CLP = "#8D6CAB", GMP = "#F4A261")

metric_label <- function(mtr) {
  if (mtr == "euclidean")
    "Euclidean distance to Heart_Granulocytes\n(lower = more similar)"
  else
    "Bray-Curtis dissimilarity to Heart_Granulocytes\n(lower = more similar)"
}

# --- Plot A: ranked, single panel ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |> dplyr::rename(value = all_of(mtr))
  src_order <- d |> dplyr::group_by(source) |>
    dplyr::summarise(m = mean(value), .groups = "drop") |>
    dplyr::arrange(m) |> dplyr::pull(source)
  
  agg <- d |> dplyr::group_by(source, population) |>
    dplyr::summarise(mean_d = mean(value),
                     sem    = sd(value) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(source = factor(source, levels = src_order))
  
  ft <- friedman_df |> dplyr::filter(metric == mtr)
  subtitle <- if (!is.na(ft$p_value))
    sprintf("Friedman across %d complete-block sources (n=%d mice): chi2 = %.2f, p = %.3g",
            ft$n_sources_tested, ft$n_mice, ft$statistic, ft$p_value)
  else "Friedman not run (fewer than 2 complete-block sources)"
  
  d_pts <- d |> dplyr::mutate(source = factor(source, levels = src_order))
  
  p <- ggplot(agg, aes(x = source, y = mean_d, fill = population)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.78) +
    geom_errorbar(aes(ymin = mean_d - sem, ymax = mean_d + sem),
                  width = 0.3, linewidth = 0.4) +
    geom_point(data = d_pts, inherit.aes = FALSE,
               mapping = aes(x = source, y = value, color = factor(mouse)),
               size = 1.9, alpha = 0.9,
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    scale_fill_manual(values = pop_palette, name = "Population") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Which (population x bone) source is closest to Heart_Granulocytes?",
         subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05a_Heart_Granulocytes_distance_ranked_", mtr),
            w = max(11, 0.55 * length(src_order) + 3), h = 5.8)
}

# --- Plot B: faceted by population, lines connecting mice across bones ---
for (mtr in c("euclidean", "bray_curtis")) {
  d <- dist_df |>
    dplyr::rename(value = all_of(mtr)) |>
    dplyr::mutate(bone = factor(bone, levels = bones))
  
  p <- d |>
    ggplot(aes(x = bone, y = value)) +
    stat_summary(aes(fill = population), fun = mean, geom = "col",
                 color = "black", linewidth = 0.3, alpha = 0.75) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.3, linewidth = 0.4) +
    geom_line(aes(group = mouse, color = factor(mouse)), alpha = 0.6, linewidth = 0.4) +
    geom_point(aes(color = factor(mouse)), size = 2, alpha = 0.9) +
    facet_wrap(~ population, nrow = 1) +
    scale_fill_manual(values = pop_palette, guide = "none") +
    scale_color_brewer(palette = "Set1", name = "Mouse") +
    labs(x = NULL, y = metric_label(mtr),
         title = "Distance to Heart_Granulocytes by bone, faceted by population") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
  
  save_plot(p, paste0("05b_Heart_Granulocytes_distance_faceted_", mtr),
            w = max(13, 3.4 * length(populations)), h = 5.2)
}

message("Done. Outputs in: ", output_dir)


