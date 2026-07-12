# combine_and_analyze_CT.R
# STEP 1: Combine C-region and T-region wide barcode tables into a single
#         CT wide table. Samples are matched between C and T by stripping
#         BOTH the _S## suffix AND the _CA / _TA region tag (since C-amplicon
#         and T-amplicon sequencing of the same physical sample carry those
#         tags). Barcodes are prefixed "C:" or "T:".
# STEP 2: Per-mouse distance-to-destination analyses on the combined matrix.
#         For each sample, C-barcodes and T-barcodes are normalized SEPARATELY
#         to sum 1 within their region, then halved and stacked so C and T
#         contribute equally to the combined frequency vector regardless of
#         their per-region sequencing depth. Samples missing one region
#         entirely get the present region renormalized to full weight.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- USER SETTINGS ------------------------------------------------------
c_file <- "C:/Users/fc809/Downloads/C_regions_results/combined_barcodes_C_regions.tsv"
t_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"

output_dir    <- "C:/Users/fc809/Downloads/CT_regions_results/"
combined_file <- file.path(output_dir, "combined_barcodes_wide_CT.tsv")

destinations <- tibble::tribble(
  ~short,                 ~pattern,
  "Thymus_DN1",           "Thymus_DN1",
  "Blood_T_cells",        "Blood_T_cells",
  "Blood_Granulocytes",   "Blood_Granulocytes",
  "Heart_Granulocytes",   "Heart_Granulocytes"
)

populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")
# -------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# STEP 1: COMBINE C AND T
# =========================================================================
read_and_prep <- function(path, region_prefix, strip_region_tag) {
  # strip_region_tag: e.g. "CA" for C-amplicon file or "TA" for T-amplicon file.
  # Removes BOTH the trailing _S## and the _<tag> suffix so the same physical
  # sample aligns between the C and T files.
  if (!file.exists(path)) stop("Input file not found: ", path)
  d <- readr::read_tsv(path, show_col_types = FALSE)
  summary_cols <- c("barcode", "n_samples_detected", "total_umi")
  sample_cols  <- setdiff(colnames(d), summary_cols)

  new_names <- sample_cols |>
    stringr::str_remove("_S\\d+$") |>
    stringr::str_remove(paste0("_", strip_region_tag, "$"))

  n_unchanged <- sum(new_names == stringr::str_remove(sample_cols, "_S\\d+$"))
  if (n_unchanged > 0) {
    warning(n_unchanged, " sample(s) in ", region_prefix,
            " file did not end with _", strip_region_tag,
            "; left as-is. Example: ",
            sample_cols[new_names == stringr::str_remove(sample_cols, "_S\\d+$")][1])
  }
  dups <- new_names[duplicated(new_names)]
  if (length(dups) > 0) {
    warning("After stripping, duplicate sample names in ", region_prefix,
            " file: ", paste(unique(dups), collapse = ", "))
  }

  d |>
    dplyr::select(barcode, dplyr::all_of(sample_cols)) |>
    dplyr::rename_with(~ new_names, dplyr::all_of(sample_cols)) |>
    dplyr::mutate(barcode = paste0(region_prefix, ":", barcode))
}

message("Reading C-region: ", c_file)
c_wide <- read_and_prep(c_file, "C", strip_region_tag = "CA")
message("Reading T-region: ", t_file)
t_wide <- read_and_prep(t_file, "T", strip_region_tag = "TA")

c_samples <- setdiff(colnames(c_wide), "barcode")
t_samples <- setdiff(colnames(t_wide), "barcode")
only_c <- setdiff(c_samples, t_samples)
only_t <- setdiff(t_samples, c_samples)
both   <- intersect(c_samples, t_samples)
message("Samples in C only: ", length(only_c),
        if (length(only_c)) paste0(" -> ", paste(only_c, collapse = ", ")) else "")
message("Samples in T only: ", length(only_t),
        if (length(only_t)) paste0(" -> ", paste(only_t, collapse = ", ")) else "")
message("Samples in both:   ", length(both))

all_samples <- union(c_samples, t_samples)
for (s in setdiff(all_samples, c_samples)) c_wide[[s]] <- 0L
for (s in setdiff(all_samples, t_samples)) t_wide[[s]] <- 0L
c_wide <- c_wide |> dplyr::select(barcode, dplyr::all_of(all_samples))
t_wide <- t_wide |> dplyr::select(barcode, dplyr::all_of(all_samples))

combined <- dplyr::bind_rows(c_wide, t_wide)

sample_cols_comb <- setdiff(colnames(combined), "barcode")
combined <- combined |>
  dplyr::mutate(
    n_samples_detected = rowSums(dplyr::across(dplyr::all_of(sample_cols_comb)) > 0),
    total_umi          = rowSums(dplyr::across(dplyr::all_of(sample_cols_comb)))
  ) |>
  dplyr::relocate(barcode, n_samples_detected, total_umi) |>
  dplyr::arrange(dplyr::desc(total_umi), barcode)

readr::write_tsv(combined, combined_file)
message("Wrote combined CT wide table: ", combined_file)
message("  Unique CT barcodes: ", nrow(combined))
message("  Samples:            ", length(sample_cols_comb))

# =========================================================================
# STEP 2: DISTANCE ANALYSIS
# =========================================================================
wide <- combined
summary_cols <- c("barcode", "n_samples_detected", "total_umi")
sample_cols  <- setdiff(colnames(wide), summary_cols)

# Sample names no longer carry _S## or _CA/_TA suffixes after combine.
meta <- tibble::tibble(sample = sample_cols) |>
  dplyr::mutate(
    parts       = stringr::str_split(sample, "_"),
    mouse       = suppressWarnings(as.integer(purrr::map_chr(parts, 1L))),
    population  = purrr::map_chr(parts, ~ if (length(.x) >= 2) .x[2] else NA_character_),
    bone        = purrr::map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
    sample_type = sample |> stringr::str_remove("^[^_]+_")
  ) |>
  dplyr::select(-parts)

src_meta_all <- meta |>
  dplyr::filter(!is.na(mouse),
                population %in% populations,
                bone        %in% bones)

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

euclidean   <- function(a, b) sqrt(sum((a - b)^2))
bray_curtis <- function(a, b) {
  denom <- sum(a + b)
  if (denom == 0) NA_real_ else sum(abs(a - b)) / denom
}

pop_palette <- c(LSK = "#2A9D8F", MPP4 = "#E9C46A",
                 CLP = "#8D6CAB", GMP = "#F4A261")

# Build a per-sample frequency matrix where C and T are normalized SEPARATELY
# per sample and then averaged 50/50 (full weight if one region is missing).
build_freq_matrix <- function(m_sub, sample_cols_local) {
  is_c <- stringr::str_starts(m_sub$barcode, "C:")
  M <- as.matrix(m_sub[, sample_cols_local])
  M_c <- M[is_c,  , drop = FALSE]
  M_t <- M[!is_c, , drop = FALSE]

  c_totals <- colSums(M_c); c_totals_safe <- ifelse(c_totals == 0, 1, c_totals)
  t_totals <- colSums(M_t); t_totals_safe <- ifelse(t_totals == 0, 1, t_totals)
  F_c <- sweep(M_c, 2, c_totals_safe, "/")
  F_t <- sweep(M_t, 2, t_totals_safe, "/")

  F_stacked <- rbind(F_c / 2, F_t / 2)
  final_totals <- colSums(F_stacked); final_totals[final_totals == 0] <- 1
  Fmat <- sweep(F_stacked, 2, final_totals, "/")
  rownames(Fmat) <- m_sub$barcode
  Fmat
}

run_destination <- function(dest_pattern, dest_short) {
  message("\n========== CT | Destination: ", dest_short,
          "  (pattern: '", dest_pattern, "') ==========")

  dest_meta <- meta |>
    dplyr::filter(!is.na(mouse),
                  stringr::str_detect(sample_type, dest_pattern))
  message("Found ", nrow(dest_meta), " destination samples across ",
          dplyr::n_distinct(dest_meta$mouse), " mice.")
  if (nrow(dest_meta) == 0) {
    warning("Skipping ", dest_short, ": no destination samples matched.")
    return(invisible(NULL))
  }
  if (anyDuplicated(dest_meta$mouse)) {
    warning("Multiple destination samples found in some mice for ", dest_short,
            "; using the first per mouse.")
    dest_meta <- dest_meta |> dplyr::group_by(mouse) |>
      dplyr::slice(1) |> dplyr::ungroup()
  }

  results <- list()
  for (m in sort(unique(dest_meta$mouse))) {
    dest_col <- dest_meta$sample[dest_meta$mouse == m]
    m_src    <- src_meta_all |> dplyr::filter(mouse == m)
    if (nrow(m_src) == 0) {
      message("  Mouse ", m, ": destination found but no qualifying sources; skipping.")
      next
    }
    m_cols <- c(dest_col, m_src$sample)

    m_sub <- wide[, c("barcode", m_cols)]
    keep  <- rowSums(as.matrix(m_sub[, m_cols])) > 0
    m_sub <- m_sub[keep, , drop = FALSE]

    Fmat <- build_freq_matrix(m_sub, m_cols)

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
  if (nrow(dist_df) == 0) {
    warning("No distances computed for ", dest_short)
    return(invisible(NULL))
  }
  dist_df <- dist_df |>
    dplyr::group_by(mouse, population, bone, source) |>
    dplyr::summarise(euclidean   = mean(euclidean),
                     bray_curtis = mean(bray_curtis),
                     .groups     = "drop")

  readr::write_csv(dist_df,
                   file.path(output_dir,
                             paste0("CT_", dest_short, "_distances_per_mouse.csv")))

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
  readr::write_csv(friedman_df,
                   file.path(output_dir,
                             paste0("CT_", dest_short, "_distance_friedman.csv")))
  print(friedman_df)

  metric_label <- function(mtr) {
    if (mtr == "euclidean")
      paste0("Euclidean distance to ", dest_short,
             "\n(C and T frequencies averaged per sample; lower = more similar)")
    else
      paste0("Bray-Curtis dissimilarity to ", dest_short,
             "\n(C and T frequencies averaged per sample; lower = more similar)")
  }

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
           title = paste0("CT barcodes: which source is closest to ", dest_short, "?"),
           subtitle = subtitle) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
            panel.grid.major.x = element_blank())

    save_plot(p, paste0("05a_CT_", dest_short, "_distance_ranked_", mtr),
              w = max(11, 0.55 * length(src_order) + 3), h = 5.8)
  }

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
      geom_line(aes(group = mouse, color = factor(mouse)),
                alpha = 0.6, linewidth = 0.4) +
      geom_point(aes(color = factor(mouse)), size = 2, alpha = 0.9) +
      facet_wrap(~ population, nrow = 1) +
      scale_fill_manual(values = pop_palette, guide = "none") +
      scale_color_brewer(palette = "Set1", name = "Mouse") +
      labs(x = NULL, y = metric_label(mtr),
           title = paste0("CT barcodes: distance to ", dest_short,
                          " by bone, faceted by population")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            strip.text = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
    save_plot(p, paste0("05b_CT_", dest_short, "_distance_faceted_", mtr),
              w = max(13, 3.4 * length(populations)), h = 5.2)
  }
}

for (i in seq_len(nrow(destinations))) {
  run_destination(destinations$pattern[i], destinations$short[i])
}

message("\nAll destinations done (CT combined). Outputs in: ", output_dir)
