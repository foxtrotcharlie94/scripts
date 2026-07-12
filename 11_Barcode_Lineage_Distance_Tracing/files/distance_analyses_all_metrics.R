# distance_analyses_all_metrics.R
# Per-mouse distance-to-destination analyses for 4 destinations and 4 metrics
# (Euclidean, Bray-Curtis, cosine distance, Jensen-Shannon distance), run for
# three datasets:
#   - C only      (C-region barcode wide table; straight per-sample normalization)
#   - T only      (T-region barcode wide table; straight per-sample normalization)
#   - CT combined (both regions; C and T normalized SEPARATELY per sample and
#                  averaged 50/50, as in combine_and_analyze_CT.R)
# Per (dataset, destination, metric): per-mouse distances CSV, Friedman CSV,
# ranked plot (05a) and faceted-by-population plot (05b), saved as PNG + PDF.
# Filenames are prefixed with C_, T_, or CT_ to keep outputs separated.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- USER SETTINGS ------------------------------------------------------
c_wide_file  <- "C:/Users/fc809/Downloads/C_regions_results/combined_barcodes_C_regions.tsv"
t_wide_file  <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"
ct_wide_file <- "C:/Users/fc809/Downloads/CT_regions_results/combined_barcodes_wide_CT.tsv"

c_output_dir  <- "C:/Users/fc809/Downloads/C_regions_results/"
t_output_dir  <- "C:/Users/fc809/Downloads/T_region_result/"
ct_output_dir <- "C:/Users/fc809/Downloads/CT_regions_results/"

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

# ---- DISTANCE FUNCTIONS -------------------------------------------------
euclidean <- function(a, b) sqrt(sum((a - b)^2))

bray_curtis <- function(a, b) {
  d <- sum(a + b)
  if (d == 0) NA_real_ else sum(abs(a - b)) / d
}

cosine_dist <- function(a, b) {
  na <- sqrt(sum(a^2)); nb <- sqrt(sum(b^2))
  if (na == 0 || nb == 0) NA_real_ else 1 - sum(a * b) / (na * nb)
}

# Jensen-Shannon distance using log2: bounded [0, 1] for two probability dists.
# Convention: 0 * log(0/x) = 0 (use only entries where both arguments are > 0).
jensen_shannon <- function(p, q) {
  m <- (p + q) / 2
  safe_kl <- function(x, y) {
    nz <- x > 0 & y > 0
    if (!any(nz)) 0 else sum(x[nz] * log2(x[nz] / y[nz]))
  }
  jsd <- 0.5 * safe_kl(p, m) + 0.5 * safe_kl(q, m)
  sqrt(max(0, jsd))
}

metric_fns <- list(
  euclidean      = euclidean,
  bray_curtis    = bray_curtis,
  cosine         = cosine_dist,
  jensen_shannon = jensen_shannon
)

metric_labels <- c(
  euclidean      = "Euclidean distance",
  bray_curtis    = "Bray-Curtis dissimilarity",
  cosine         = "Cosine distance (1 - similarity)",
  jensen_shannon = "Jensen-Shannon distance"
)

pop_palette <- c(LSK = "#2A9D8F", MPP4 = "#E9C46A",
                 CLP = "#8D6CAB", GMP = "#F4A261")

save_plot <- function(plot, name, output_dir, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# ---- MAIN DATASET RUNNER ------------------------------------------------
run_dataset <- function(wide_file, output_dir, dataset_name, use_ct_norm) {
  message("\n###############################################")
  message("# Dataset: ", dataset_name)
  message("# Input:   ", wide_file)
  message("# Output:  ", output_dir)
  message("###############################################")

  if (!file.exists(wide_file)) {
    warning("Input file not found; skipping ", dataset_name, ": ", wide_file)
    return(invisible(NULL))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  wide <- readr::read_tsv(wide_file, show_col_types = FALSE)
  summary_cols <- c("barcode", "n_samples_detected", "total_umi")
  sample_cols  <- setdiff(colnames(wide), summary_cols)

  meta <- tibble::tibble(sample = sample_cols) |>
    dplyr::mutate(
      parts       = stringr::str_split(sample, "_"),
      mouse       = suppressWarnings(as.integer(purrr::map_chr(parts, 1L))),
      population  = purrr::map_chr(parts, ~ if (length(.x) >= 2) .x[2] else NA_character_),
      bone        = purrr::map_chr(parts, ~ if (length(.x) >= 3) .x[3] else NA_character_),
      sample_type = sample |>
        stringr::str_remove("^[^_]+_") |>
        stringr::str_remove("_S\\d+$")
    ) |>
    dplyr::select(-parts)

  src_meta_all <- meta |>
    dplyr::filter(!is.na(mouse),
                  population %in% populations,
                  bone        %in% bones)

  # Frequency-matrix builder for a per-mouse subset
  build_freq_for_subset <- function(m_sub, sample_cols_local) {
    if (use_ct_norm) {
      is_c_sub <- stringr::str_starts(m_sub$barcode, "C:")
      M <- as.matrix(m_sub[, sample_cols_local])
      M_c <- M[is_c_sub,  , drop = FALSE]
      M_t <- M[!is_c_sub, , drop = FALSE]
      c_totals <- colSums(M_c); c_totals[c_totals == 0] <- 1
      t_totals <- colSums(M_t); t_totals[t_totals == 0] <- 1
      F_c <- sweep(M_c, 2, c_totals, "/")
      F_t <- sweep(M_t, 2, t_totals, "/")
      F_stacked <- rbind(F_c / 2, F_t / 2)
      totals <- colSums(F_stacked); totals[totals == 0] <- 1
      Fmat <- sweep(F_stacked, 2, totals, "/")
      bc_order <- c(m_sub$barcode[is_c_sub], m_sub$barcode[!is_c_sub])
      rownames(Fmat) <- bc_order
      Fmat[m_sub$barcode, , drop = FALSE]
    } else {
      M <- as.matrix(m_sub[, sample_cols_local])
      totals <- colSums(M); totals[totals == 0] <- 1
      Fmat <- sweep(M, 2, totals, "/")
      rownames(Fmat) <- m_sub$barcode
      Fmat
    }
  }

  # Per-destination work
  for (di in seq_len(nrow(destinations))) {
    dest_pattern <- destinations$pattern[di]
    dest_short   <- destinations$short[di]
    message("\n========== ", dataset_name, " | Destination: ", dest_short, " ==========")

    dest_meta <- meta |>
      dplyr::filter(!is.na(mouse),
                    stringr::str_detect(sample_type, dest_pattern))
    if (nrow(dest_meta) == 0) {
      warning("No destination samples matched '", dest_pattern, "' in ", dataset_name)
      next
    }
    if (anyDuplicated(dest_meta$mouse)) {
      warning("Multiple destination samples per mouse for ", dest_short,
              " in ", dataset_name, "; using the first per mouse.")
      dest_meta <- dest_meta |> dplyr::group_by(mouse) |>
        dplyr::slice(1) |> dplyr::ungroup()
    }

    # Per-mouse distance computation
    results <- list()
    for (m in sort(unique(dest_meta$mouse))) {
      dest_col <- dest_meta$sample[dest_meta$mouse == m]
      m_src    <- src_meta_all |> dplyr::filter(mouse == m)
      if (nrow(m_src) == 0) next
      m_cols <- c(dest_col, m_src$sample)

      m_sub <- wide[, c("barcode", m_cols)]
      keep  <- rowSums(as.matrix(m_sub[, m_cols])) > 0
      m_sub <- m_sub[keep, , drop = FALSE]

      Fmat <- build_freq_for_subset(m_sub, m_cols)
      dest_freq <- Fmat[, dest_col]

      for (i in seq_len(nrow(m_src))) {
        src_col  <- m_src$sample[i]
        src_freq <- Fmat[, src_col]
        row <- tibble::tibble(
          mouse      = m,
          population = m_src$population[i],
          bone       = m_src$bone[i],
          source     = paste(m_src$population[i], m_src$bone[i], sep = "_"),
          sample     = src_col
        )
        for (mtr in names(metric_fns)) {
          row[[mtr]] <- metric_fns[[mtr]](src_freq, dest_freq)
        }
        results[[length(results) + 1]] <- row
      }
    }
    dist_df <- dplyr::bind_rows(results)
    if (nrow(dist_df) == 0) {
      warning("No distances computed for ", dest_short, " in ", dataset_name)
      next
    }
    dist_df <- dist_df |>
      dplyr::group_by(mouse, population, bone, source) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(names(metric_fns)), mean),
                       .groups = "drop")

    readr::write_csv(dist_df,
      file.path(output_dir,
                paste0(dataset_name, "_", dest_short, "_distances_per_mouse.csv")))

    # Friedman per metric
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
    friedman_df <- dplyr::bind_rows(lapply(names(metric_fns), friedman_one))
    readr::write_csv(friedman_df,
      file.path(output_dir,
                paste0(dataset_name, "_", dest_short, "_distance_friedman.csv")))
    print(friedman_df)

    # Plots per metric
    for (mtr in names(metric_fns)) {
      ylab_text <- paste0(metric_labels[[mtr]], " to ", dest_short,
                          "\n(lower = more similar)")
      d <- dist_df |> dplyr::rename(value = all_of(mtr))

      # 05a ranked
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
        labs(x = NULL, y = ylab_text,
             title = paste0(dataset_name, " barcodes: which source is closest to ",
                            dest_short, "?"),
             subtitle = subtitle) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
              panel.grid.major.x = element_blank())
      save_plot(p,
        paste0("05a_", dataset_name, "_", dest_short, "_distance_ranked_", mtr),
        output_dir,
        w = max(11, 0.55 * length(src_order) + 3), h = 5.8)

      # 05b faceted by population
      d_facet <- d |> dplyr::mutate(bone = factor(bone, levels = bones))
      p <- d_facet |>
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
        labs(x = NULL, y = ylab_text,
             title = paste0(dataset_name, " barcodes: distance to ", dest_short,
                            " by bone, faceted by population")) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1),
              strip.text = element_text(face = "bold"),
              panel.grid.major.x = element_blank())
      save_plot(p,
        paste0("05b_", dataset_name, "_", dest_short, "_distance_faceted_", mtr),
        output_dir,
        w = max(13, 3.4 * length(populations)), h = 5.2)
    }
  }
}

# ---- RUN ALL THREE DATASETS --------------------------------------------
run_dataset(c_wide_file,  c_output_dir,  "C",  use_ct_norm = FALSE)
run_dataset(t_wide_file,  t_output_dir,  "T",  use_ct_norm = FALSE)
run_dataset(ct_wide_file, ct_output_dir, "CT", use_ct_norm = TRUE)

message("\nAll datasets done.")
