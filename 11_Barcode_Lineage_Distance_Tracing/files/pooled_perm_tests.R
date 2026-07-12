# pooled_perm_tests.R
# Pool barcodes across mice per sample_type, then run UMI-level Binomial
# permutation tests for source preference between user-specified destination
# pairs. Now runs FOUR distance metrics per contrast:
#   - Euclidean
#   - Bray-Curtis
#   - cosine distance (1 - cosine similarity)
#   - Jensen-Shannon distance (log2)
#
# Per (contrast, metric) outputs:
#   CT_pooled_<pos>_vs_<neg>_perm_<metric>.csv
#   07_CT_pooled_<pos>_vs_<neg>_perm_<metric>.{png,pdf}

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- USER SETTINGS ------------------------------------------------------
combined_file <- "C:/Users/fc809/Downloads/CT_regions_results/combined_barcodes_wide_CT.tsv"
output_dir    <- "C:/Users/fc809/Downloads/CT_regions_results/"

# Each row defines one contrast. Positive delta = closer to dest_pos.
contrasts <- tibble::tribble(
  ~dest_pos_pattern,   ~dest_neg_pattern,
  "Blood_T_cells",     "Blood_Granulocytes",
  "Thymus_DN1",        "Blood_Granulocytes"
  # add more rows here for additional contrasts
)

populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")

n_perms <- 5000
seed    <- 1
# -------------------------------------------------------------------------

set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(combined_file)) {
  stop("Combined CT file not found: ", combined_file,
       "\nRun combine_and_analyze_CT.R first.")
}

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
  euclidean      = "Euclidean",
  bray_curtis    = "Bray-Curtis",
  cosine         = "cosine",
  jensen_shannon = "Jensen-Shannon"
)

# ---- LOAD & POOL --------------------------------------------------------
message("Reading: ", combined_file)
wide <- readr::read_tsv(combined_file, show_col_types = FALSE)
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

long <- wide |>
  dplyr::select(barcode, dplyr::all_of(sample_cols)) |>
  tidyr::pivot_longer(-barcode, names_to = "sample", values_to = "umi") |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::filter(!is.na(mouse))

pooled_long <- long |>
  dplyr::group_by(barcode, sample_type, population, bone) |>
  dplyr::summarise(umi = sum(umi), .groups = "drop")

# ---- HELPERS ------------------------------------------------------------
save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

# Build per-sample frequency matrix with C/T-region averaging
make_freq_builder <- function(pooled_wide_local) {
  is_c <- stringr::str_starts(pooled_wide_local$barcode, "C:")
  bc_order <- pooled_wide_local$barcode
  bc_split_order <- c(bc_order[is_c], bc_order[!is_c])

  function(counts_mat) {
    M_c <- counts_mat[is_c,  , drop = FALSE]
    M_t <- counts_mat[!is_c, , drop = FALSE]
    c_totals <- colSums(M_c); c_totals[c_totals == 0] <- 1
    t_totals <- colSums(M_t); t_totals[t_totals == 0] <- 1
    F_c <- sweep(M_c, 2, c_totals, "/")
    F_t <- sweep(M_t, 2, t_totals, "/")
    F_stacked <- rbind(F_c / 2, F_t / 2)
    totals <- colSums(F_stacked); totals[totals == 0] <- 1
    Fmat <- sweep(F_stacked, 2, totals, "/")
    rownames(Fmat) <- bc_split_order
    Fmat[bc_order, , drop = FALSE]
  }
}

# ---- RUN ONE CONTRAST ---------------------------------------------------
run_contrast <- function(dest_pos_pattern, dest_neg_pattern) {
  message("\n#####################################################")
  message("# Contrast: ", dest_pos_pattern, "  (+)  vs  ", dest_neg_pattern, "  (-)")
  message("#####################################################")

  dest_pos_types <- unique(pooled_long$sample_type[
    stringr::str_detect(pooled_long$sample_type, dest_pos_pattern)])
  dest_neg_types <- unique(pooled_long$sample_type[
    stringr::str_detect(pooled_long$sample_type, dest_neg_pattern)])
  if (length(dest_pos_types) != 1) {
    stop("Expected exactly one sample_type matching '", dest_pos_pattern,
         "'; got ", length(dest_pos_types), ": ",
         paste(dest_pos_types, collapse = ", "))
  }
  if (length(dest_neg_types) != 1) {
    stop("Expected exactly one sample_type matching '", dest_neg_pattern,
         "'; got ", length(dest_neg_types), ": ",
         paste(dest_neg_types, collapse = ", "))
  }
  message("  dest_pos: ", dest_pos_types)
  message("  dest_neg: ", dest_neg_types)

  src_types <- pooled_long |>
    dplyr::filter(population %in% populations, bone %in% bones) |>
    dplyr::distinct(sample_type, population, bone) |>
    dplyr::arrange(sample_type)
  message("  Source sample_types: ", nrow(src_types))

  pooled_wide <- pooled_long |>
    dplyr::select(-population, -bone) |>
    tidyr::pivot_wider(names_from = sample_type, values_from = umi, values_fill = 0)

  relevant_cols <- c(dest_pos_types, dest_neg_types, src_types$sample_type)
  relevant_cols <- intersect(relevant_cols, colnames(pooled_wide))
  mask <- rowSums(as.matrix(pooled_wide[, relevant_cols])) > 0
  pooled_wide <- pooled_wide[mask, , drop = FALSE]
  message("  Barcodes after filtering for relevance: ", nrow(pooled_wide))

  build_freq <- make_freq_builder(pooled_wide)

  # Observed
  M_obs <- as.matrix(pooled_wide[, relevant_cols])
  F_obs <- build_freq(M_obs)
  T_freq_obs <- F_obs[, dest_pos_types]
  G_freq_obs <- F_obs[, dest_neg_types]
  src_freqs  <- F_obs[, src_types$sample_type, drop = FALSE]
  n_src <- ncol(src_freqs)

  # Observed deltas per metric
  obs_deltas <- list()
  for (mtr in names(metric_fns)) {
    fn <- metric_fns[[mtr]]
    d <- apply(src_freqs, 2, function(s) c(d_T = fn(s, T_freq_obs),
                                           d_G = fn(s, G_freq_obs)))
    obs_deltas[[mtr]] <- d["d_G", ] - d["d_T", ]
  }

  # Permutation
  counts_T <- pooled_wide[[dest_pos_types]]
  counts_G <- pooled_wide[[dest_neg_types]]
  N_total <- counts_T + counts_G
  N_T <- sum(counts_T); N_G <- sum(counts_G)
  p_T <- N_T / (N_T + N_G)
  message(sprintf("  Total T UMIs: %d  |  Total G UMIs: %d  |  p_T = %.4f",
                  N_T, N_G, p_T))

  null_deltas <- lapply(metric_fns, function(fn) {
    matrix(NA_real_, n_perms, n_src,
           dimnames = list(NULL, colnames(src_freqs)))
  })

  message("  Running ", n_perms, " permutations across ",
          length(metric_fns), " metrics...")
  pb <- txtProgressBar(min = 0, max = n_perms, style = 3)
  for (i in seq_len(n_perms)) {
    T_p <- rbinom(length(N_total), N_total, p_T)
    G_p <- N_total - T_p
    M_p <- cbind(T = T_p, G = G_p)
    Fp  <- build_freq(M_p)
    T_F <- Fp[, "T"]; G_F <- Fp[, "G"]
    for (mtr in names(metric_fns)) {
      fn <- metric_fns[[mtr]]
      for (j in seq_len(n_src)) {
        s <- src_freqs[, j]
        null_deltas[[mtr]][i, j] <- fn(s, G_F) - fn(s, T_F)
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Build results + plots per metric
  for (mtr in names(metric_fns)) {
    obs <- obs_deltas[[mtr]]
    nullmat <- null_deltas[[mtr]]
    p_two <- vapply(seq_along(obs),
                    function(j) mean(abs(nullmat[, j]) >= abs(obs[j]), na.rm = TRUE),
                    numeric(1))
    null_lo <- apply(nullmat, 2, quantile, probs = 0.025, na.rm = TRUE)
    null_hi <- apply(nullmat, 2, quantile, probs = 0.975, na.rm = TRUE)
    res <- tibble::tibble(
      sample_type = src_types$sample_type,
      population  = src_types$population,
      bone        = src_types$bone,
      mean_delta  = obs,
      p_two_sided = p_two,
      p_bh        = p.adjust(p_two, method = "BH"),
      null_p2_5   = null_lo,
      null_p97_5  = null_hi
    ) |>
      dplyr::arrange(mean_delta)

    csv_name <- paste0("CT_pooled_", dest_pos_pattern, "_vs_",
                       dest_neg_pattern, "_perm_", mtr, ".csv")
    readr::write_csv(res, file.path(output_dir, csv_name))

    # Plot
    pos_label <- paste0("closer to ", dest_pos_pattern)
    neg_label <- paste0("closer to ", dest_neg_pattern)
    fill_vals <- setNames(c("#2A9D8F", "#E76F51"), c(pos_label, neg_label))

    res_plot <- res |>
      dplyr::mutate(
        direction = ifelse(mean_delta > 0, pos_label, neg_label),
        stars = dplyr::case_when(
          p_bh < 0.001 ~ "***",
          p_bh < 0.01  ~ "**",
          p_bh < 0.05  ~ "*",
          TRUE ~ ""
        ),
        label = ifelse(stars == "",
                       sprintf("p=%.3g", p_two_sided),
                       sprintf("%s p=%.3g (BH %.3g)", stars, p_two_sided, p_bh))
      )

    ylab_str <- sprintf("Mean %s delta = d(. -> %s) - d(. -> %s)   [pooled across mice]",
                        metric_labels[[mtr]], dest_neg_pattern, dest_pos_pattern)

    p <- ggplot(res_plot, aes(x = reorder(sample_type, mean_delta),
                              y = mean_delta, fill = direction)) +
      geom_col(color = "black", linewidth = 0.3, alpha = 0.82) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      geom_errorbar(aes(ymin = null_p2_5, ymax = null_p97_5),
                    width = 0.25, linewidth = 0.3, color = "grey30") +
      geom_text(aes(label = label,
                    hjust = ifelse(mean_delta > 0, -0.05, 1.05)),
                size = 3) +
      scale_fill_manual(values = fill_vals, name = NULL) +
      coord_flip() +
      labs(x = NULL, y = ylab_str,
           title = sprintf("Pooled %s vs %s (metric: %s)",
                           dest_pos_pattern, dest_neg_pattern, metric_labels[[mtr]]),
           subtitle = sprintf("UMI-level Binomial permutation (n_perms = %d). Grey bars = 95%% null interval.",
                              n_perms)) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.major.y = element_blank(),
            plot.subtitle = element_text(size = 9))

    plot_name <- sprintf("07_CT_pooled_%s_vs_%s_perm_%s",
                         dest_pos_pattern, dest_neg_pattern, mtr)
    save_plot(p, plot_name, w = 13, h = max(5, 0.32 * nrow(res) + 2))
  }

  message("  Contrast done: ", dest_pos_pattern, " vs ", dest_neg_pattern)
}

# ---- RUN ALL CONTRASTS --------------------------------------------------
for (k in seq_len(nrow(contrasts))) {
  run_contrast(contrasts$dest_pos_pattern[k], contrasts$dest_neg_pattern[k])
}

message("\nAll contrasts done. Outputs in: ", output_dir)
