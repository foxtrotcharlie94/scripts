# pooled_perm_test_TvG.R
# Self-contained script:
#   STEP A: read C-region and T-region wide barcode tables, combine into one CT
#           wide table (samples matched by stripping _S## suffix; barcodes
#           prefixed "C:" / "T:"), save to output_dir.
#   STEP B: pool UMIs across mice per sample_type and run a UMI-level
#           Binomial permutation test for source preference between two
#           destinations (default: Blood_T_cells vs Blood_Granulocytes).

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- USER SETTINGS ------------------------------------------------------
# Inputs: the two per-region combined wide tables.
c_file <- "C:/Users/fc809/Downloads/C_regions_results/combined_barcodes_C_regions.tsv"
t_file <- "C:/Users/fc809/Downloads/T_region_result/combined_barcodes_T_region.tsv"

output_dir <- "C:/Users/fc809/Downloads/CT_regions_results/"
combined_file <- file.path(output_dir, "combined_barcodes_wide_CT.tsv")

dest_pos_pattern <- "Blood_T_cells"          # positive delta direction
dest_neg_pattern <- "Blood_Granulocytes"     # negative delta direction

populations <- c("LSK", "MPP4", "CLP", "GMP")
bones       <- c("Femur", "Spine", "Hip", "Tibia", "Skull", "Sternum", "Pelvis")

n_perms <- 5000
seed    <- 1
# -------------------------------------------------------------------------

set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# STEP A: COMBINE C AND T (writes combined_barcodes_wide_CT.tsv to output_dir)
# =========================================================================
read_and_prep <- function(path, region_prefix, strip_region_tag) {
  # strip_region_tag: e.g. "CA" for C-amplicon file or "TA" for T-amplicon file.
  # Removes BOTH the trailing _S## and the _<tag> suffix so that the same
  # physical sample aligns between the C and T files.
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

message("STEP A: combining C and T region wide tables.")
message("  Reading C-region: ", c_file)
c_wide <- read_and_prep(c_file, "C", strip_region_tag = "CA")
message("  Reading T-region: ", t_file)
t_wide <- read_and_prep(t_file, "T", strip_region_tag = "TA")

c_samples <- setdiff(colnames(c_wide), "barcode")
t_samples <- setdiff(colnames(t_wide), "barcode")
only_c <- setdiff(c_samples, t_samples)
only_t <- setdiff(t_samples, c_samples)
both   <- intersect(c_samples, t_samples)
message("  Samples in C only: ", length(only_c))
message("  Samples in T only: ", length(only_t))
message("  Samples in both:   ", length(both))

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
message("  Wrote: ", combined_file)
message("  Unique CT barcodes: ", nrow(combined))
message("  Samples:            ", length(sample_cols_comb))

# =========================================================================
# STEP B: POOLED PERMUTATION TEST
# =========================================================================
wide <- combined
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

# Pool UMIs across mice per (barcode, sample_type)
long <- wide |>
  dplyr::select(barcode, dplyr::all_of(sample_cols)) |>
  tidyr::pivot_longer(-barcode, names_to = "sample", values_to = "umi") |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::filter(!is.na(mouse))

pooled_long <- long |>
  dplyr::group_by(barcode, sample_type, population, bone) |>
  dplyr::summarise(umi = sum(umi), .groups = "drop")

dest_pos_types <- unique(pooled_long$sample_type[
  stringr::str_detect(pooled_long$sample_type, dest_pos_pattern)])
dest_neg_types <- unique(pooled_long$sample_type[
  stringr::str_detect(pooled_long$sample_type, dest_neg_pattern)])
stopifnot("expected exactly one matching sample_type for dest_pos" =
            length(dest_pos_types) == 1,
          "expected exactly one matching sample_type for dest_neg" =
            length(dest_neg_types) == 1)
message("\nSTEP B: pooled permutation test")
message("  dest_pos sample_type: ", dest_pos_types)
message("  dest_neg sample_type: ", dest_neg_types)

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

is_c <- stringr::str_starts(pooled_wide$barcode, "C:")

build_freq <- function(counts_mat) {
  M_c <- counts_mat[is_c,  , drop = FALSE]
  M_t <- counts_mat[!is_c, , drop = FALSE]
  c_totals <- colSums(M_c); c_totals[c_totals == 0] <- 1
  t_totals <- colSums(M_t); t_totals[t_totals == 0] <- 1
  F_c <- sweep(M_c, 2, c_totals, "/")
  F_t <- sweep(M_t, 2, t_totals, "/")
  F_stacked <- rbind(F_c / 2, F_t / 2)
  totals <- colSums(F_stacked); totals[totals == 0] <- 1
  Fmat <- sweep(F_stacked, 2, totals, "/")
  bc_order <- c(pooled_wide$barcode[is_c], pooled_wide$barcode[!is_c])
  rownames(Fmat) <- bc_order
  Fmat[pooled_wide$barcode, , drop = FALSE]
}

M_obs <- as.matrix(pooled_wide[, relevant_cols])
F_obs <- build_freq(M_obs)
T_freq_obs <- F_obs[, dest_pos_types]
G_freq_obs <- F_obs[, dest_neg_types]
src_freqs  <- F_obs[, src_types$sample_type, drop = FALSE]

euclidean   <- function(a, b) sqrt(sum((a - b)^2))
bray_curtis <- function(a, b) {
  d <- sum(a + b)
  if (d == 0) NA_real_ else sum(abs(a - b)) / d
}
dist_to_TG <- function(src_F, T_F, G_F, fn) {
  apply(src_F, 2, function(s) c(d_T = fn(s, T_F), d_G = fn(s, G_F)))
}
obs_dists_euc <- dist_to_TG(src_freqs, T_freq_obs, G_freq_obs, euclidean)
obs_dists_bc  <- dist_to_TG(src_freqs, T_freq_obs, G_freq_obs, bray_curtis)
obs_delta_euc <- obs_dists_euc["d_G", ] - obs_dists_euc["d_T", ]
obs_delta_bc  <- obs_dists_bc["d_G", ]  - obs_dists_bc["d_T", ]

counts_T <- pooled_wide[[dest_pos_types]]
counts_G <- pooled_wide[[dest_neg_types]]
N_total  <- counts_T + counts_G
N_T <- sum(counts_T); N_G <- sum(counts_G)
p_T <- N_T / (N_T + N_G)
message(sprintf("  Total T UMIs: %d  |  Total G UMIs: %d  |  p_T = %.4f",
                N_T, N_G, p_T))

normalize_TG <- function(T_counts, G_counts) {
  M <- cbind(T = T_counts, G = G_counts)
  M_c <- M[is_c,  , drop = FALSE]
  M_t <- M[!is_c, , drop = FALSE]
  c_totals <- colSums(M_c); c_totals[c_totals == 0] <- 1
  t_totals <- colSums(M_t); t_totals[t_totals == 0] <- 1
  F_c <- sweep(M_c, 2, c_totals, "/")
  F_t <- sweep(M_t, 2, t_totals, "/")
  F_stacked <- rbind(F_c / 2, F_t / 2)
  totals <- colSums(F_stacked); totals[totals == 0] <- 1
  Fmat <- sweep(F_stacked, 2, totals, "/")
  bc_order <- c(pooled_wide$barcode[is_c], pooled_wide$barcode[!is_c])
  rownames(Fmat) <- bc_order
  Fmat[pooled_wide$barcode, , drop = FALSE]
}

# Sanity: observed normalization matches
chk <- normalize_TG(counts_T, counts_G)
stopifnot(all.equal(unname(chk[, "T"]), unname(T_freq_obs), tolerance = 1e-10),
          all.equal(unname(chk[, "G"]), unname(G_freq_obs), tolerance = 1e-10))

n_src <- ncol(src_freqs)
null_delta_euc <- matrix(NA_real_, n_perms, n_src,
                         dimnames = list(NULL, colnames(src_freqs)))
null_delta_bc  <- matrix(NA_real_, n_perms, n_src,
                         dimnames = list(NULL, colnames(src_freqs)))

message("  Running ", n_perms, " permutations (will take a minute)...")
pb <- txtProgressBar(min = 0, max = n_perms, style = 3)
for (i in seq_len(n_perms)) {
  T_p <- rbinom(length(N_total), N_total, p_T)
  G_p <- N_total - T_p
  Fp  <- normalize_TG(T_p, G_p)
  T_F <- Fp[, "T"]; G_F <- Fp[, "G"]
  for (j in seq_len(n_src)) {
    s <- src_freqs[, j]
    null_delta_euc[i, j] <- euclidean(s, G_F)   - euclidean(s, T_F)
    null_delta_bc[i, j]  <- bray_curtis(s, G_F) - bray_curtis(s, T_F)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

build_results <- function(obs, nullmat) {
  p_two <- vapply(seq_along(obs),
                  function(j) mean(abs(nullmat[, j]) >= abs(obs[j])),
                  numeric(1))
  null_lo <- apply(nullmat, 2, quantile, probs = 0.025, na.rm = TRUE)
  null_hi <- apply(nullmat, 2, quantile, probs = 0.975, na.rm = TRUE)
  tibble::tibble(
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
}

res_euc <- build_results(obs_delta_euc, null_delta_euc)
res_bc  <- build_results(obs_delta_bc,  null_delta_bc)

readr::write_csv(res_euc,
                 file.path(output_dir, paste0("CT_pooled_", dest_pos_pattern, "_vs_",
                                              dest_neg_pattern, "_perm_euclidean.csv")))
readr::write_csv(res_bc,
                 file.path(output_dir, paste0("CT_pooled_", dest_pos_pattern, "_vs_",
                                              dest_neg_pattern, "_perm_bray_curtis.csv")))

message("\nEuclidean results:")
print(res_euc, n = nrow(res_euc))

# ---- PLOTS -------------------------------------------------------------
save_plot <- function(plot, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, bg = "white")
}

plot_tornado <- function(res, metric_name) {
  pos_label <- paste0("closer to ", dest_pos_pattern)
  neg_label <- paste0("closer to ", dest_neg_pattern)
  fill_vals <- setNames(c("#2A9D8F", "#E76F51"), c(pos_label, neg_label))
  
  res <- res |>
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
                      ifelse(metric_name == "euclidean", "Euclidean", "Bray-Curtis"),
                      dest_neg_pattern, dest_pos_pattern)
  
  ggplot(res, aes(x = reorder(sample_type, mean_delta), y = mean_delta,
                  fill = direction)) +
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
         title = sprintf("Pooled %s vs %s: per-source preference",
                         dest_pos_pattern, dest_neg_pattern),
         subtitle = sprintf("UMI-level Binomial permutation (n_perms = %d). Grey bars = 95%% null interval. BH-adjusted p shown when < 0.05.",
                            n_perms)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          plot.subtitle = element_text(size = 9))
}

p_euc <- plot_tornado(res_euc, "euclidean")
save_plot(p_euc,
          sprintf("07_CT_pooled_%s_vs_%s_perm_euclidean",
                  dest_pos_pattern, dest_neg_pattern),
          w = 13, h = max(5, 0.32 * nrow(res_euc) + 2))

p_bc <- plot_tornado(res_bc, "bray_curtis")
save_plot(p_bc,
          sprintf("07_CT_pooled_%s_vs_%s_perm_bray_curtis",
                  dest_pos_pattern, dest_neg_pattern),
          w = 13, h = max(5, 0.32 * nrow(res_bc) + 2))

message("\nDone. Outputs in: ", output_dir)