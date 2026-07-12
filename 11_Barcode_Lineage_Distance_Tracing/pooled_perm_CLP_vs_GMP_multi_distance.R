# =============================================================================
# Pooled Permutation Test — CLP versus GMP contribution to T lineage
# Metrics: Bray-Curtis, Cosine, and Jensen-Shannon
#
# Tests whether Blood T cells or Thymus DN1 are closer to CLP than GMP,
# pooled across mice, for each bone.
#
# Uses UNFILTERED barcodes, source-side bone >= 2 UMIs.
# UMIs are pooled across mice per sample type before permutation.
#
# Interpretation of delta values:
#   delta = distance(GMP, destination) - distance(CLP, destination)
#   positive delta = destination is closer to CLP than GMP
#   negative delta = destination is closer to GMP than CLP
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

base_dir      <- "/Users/karingustafsson/Documents/Jobb/2026/Results/Sternum/DARLIN"
output_dir    <- file.path(base_dir, "pooled_perm_CLP_vs_GMP_multi_distance")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

bones_of_interest <- c("Femur", "Sternum", "Pelvis", "Spine")
source_populations <- c("CLP", "GMP")
bone_min_umi      <- 2
n_perms           <- 5000
seed              <- 1
pseudo_count      <- 1e-12

# Destination populations to test.
# Each destination is compared to CLP versus GMP, not to Blood Granulocytes.
destinations <- tribble(
  ~dest_pattern,     ~dest_label,
  "Blood_T_cells",   "Blood T cells",
  "Thymus_DN1",      "DN1 Thymocytes"
)

# -----------------------------------------------------------------------------
# LOAD UNFILTERED BARCODE COUNTS
# -----------------------------------------------------------------------------

cat("Reading unfiltered C and T tables...\n")
c_wide_raw <- read_tsv(
  file.path(base_dir, "combined_barcodes_C_regions.tsv"),
  show_col_types = FALSE)
t_wide_raw <- read_tsv(
  file.path(base_dir, "combined_barcodes_T_regions.tsv"),
  show_col_types = FALSE)

strip_suffix <- function(x) {
  x <- gsub("_S\\d+$", "", x)
  x <- gsub("_(CA|TA)$", "", x)
  x
}

meta_cols <- c("barcode", "n_samples_detected", "total_umi")

c_counts <- c_wide_raw %>%
  select(barcode, all_of(setdiff(colnames(c_wide_raw), meta_cols))) %>%
  rename_with(strip_suffix, -barcode) %>%
  mutate(barcode = paste0("C:", barcode))

t_counts <- t_wide_raw %>%
  select(barcode, all_of(setdiff(colnames(t_wide_raw), meta_cols))) %>%
  rename_with(strip_suffix, -barcode) %>%
  mutate(barcode = paste0("T:", barcode))

all_samples <- union(setdiff(colnames(c_counts), "barcode"),
                     setdiff(colnames(t_counts), "barcode"))

add_missing <- function(df, all_cols) {
  missing <- setdiff(all_cols, colnames(df)[-1])
  for (col in missing) df[[col]] <- 0
  df %>% select(barcode, all_of(all_cols))
}

raw_counts <- bind_rows(
  add_missing(c_counts, all_samples),
  add_missing(t_counts, all_samples)
)

cat("Wide table:", nrow(raw_counts), "barcodes x", ncol(raw_counts), "columns\n")

# -----------------------------------------------------------------------------
# PARSE SAMPLE METADATA
# -----------------------------------------------------------------------------

sample_meta <- tibble(sample = all_samples) %>%
  mutate(
    mouse      = map_chr(sample, ~ str_split(.x, "_")[[1]][1]),
    population = map_chr(sample, ~ str_split(.x, "_")[[1]][2]),
    bone       = map_chr(sample, ~ str_split(.x, "_")[[1]][3])
  ) %>%
  filter(!is.na(suppressWarnings(as.integer(mouse)))) %>%
  mutate(bone = str_replace(bone, "Pelivis", "Pelvis"))

mice <- unique(sample_meta$mouse)
cat("Mice found:", paste(mice, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# DISTANCE FUNCTIONS
# -----------------------------------------------------------------------------

normalize <- function(v) {
  s <- sum(v)
  if (s == 0) return(v)
  v / s
}

bray_curtis_dist <- function(a, b) {
  denom <- sum(a + b)
  if (denom == 0) return(NA_real_)
  sum(abs(a - b)) / denom
}

cosine_dist <- function(a, b) {
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (denom == 0) return(NA_real_)
  1 - sum(a * b) / denom
}

js_dist <- function(p, q) {
  s_p <- sum(p)
  s_q <- sum(q)
  if (s_p == 0 || s_q == 0) return(NA_real_)
  p <- p / s_p
  q <- q / s_q
  m <- (p + q) / 2
  kl <- function(x, y) {
    idx <- x > 0 & y > 0
    sum(x[idx] * log2(x[idx] / y[idx]))
  }
  sqrt(0.5 * kl(p, m) + 0.5 * kl(q, m))
}

metric_functions <- list(
  bray_curtis    = bray_curtis_dist,
  cosine         = cosine_dist,
  jensen_shannon = js_dist
)

calc_distances <- function(source_freq, dest_freq) {
  tibble(
    metric = names(metric_functions),
    distance = map_dbl(metric_functions, ~ .x(source_freq, dest_freq))
  )
}

# -----------------------------------------------------------------------------
# POOL UMIs ACROSS MICE
# Apply bone_min_umi threshold to CLP and GMP source samples before pooling.
# -----------------------------------------------------------------------------

pool_source <- function(population_name, bone_name) {
  source_samples <- sample_meta %>%
    filter(population == population_name, bone == bone_name) %>%
    pull(sample)

  if (length(source_samples) == 0) return(rep(0, nrow(raw_counts)))

  mat <- raw_counts %>%
    select(all_of(source_samples)) %>%
    as.matrix()

  mat[mat < bone_min_umi] <- 0
  rowSums(mat)
}

pool_dest <- function(dest_pattern) {
  dest_samples <- sample_meta %>%
    filter(str_detect(sample, dest_pattern)) %>%
    pull(sample)

  if (length(dest_samples) == 0) return(rep(0, nrow(raw_counts)))

  raw_counts %>%
    select(all_of(dest_samples)) %>%
    as.matrix() %>%
    rowSums()
}

# -----------------------------------------------------------------------------
# MAIN PERMUTATION TEST
# -----------------------------------------------------------------------------

set.seed(seed)
all_results <- list()

for (di in seq_len(nrow(destinations))) {
  dest_pat <- destinations$dest_pattern[di]
  dest_lab <- destinations$dest_label[di]

  cat("\n=== Destination:", dest_lab, "===\n")

  dest_vec <- pool_dest(dest_pat)
  dest_freq <- normalize(dest_vec)

  cat("  Destination barcodes:", sum(dest_freq > 0), "\n")

  if (sum(dest_freq) == 0) {
    cat("  No destination UMIs found -- skipping\n")
    next
  }

  for (bone_name in bones_of_interest) {
    cat("  Bone:", bone_name, "\n")

    clp_vec <- pool_source("CLP", bone_name)
    gmp_vec <- pool_source("GMP", bone_name)

    n_clp_used <- sum(clp_vec > 0)
    n_gmp_used <- sum(gmp_vec > 0)

    if (n_clp_used == 0 || n_gmp_used == 0) {
      cat("    Missing CLP or GMP barcodes >=", bone_min_umi, "UMIs -- skipping\n")
      next
    }

    clp_freq <- normalize(clp_vec)
    gmp_freq <- normalize(gmp_vec)

    clp_dist <- calc_distances(clp_freq, dest_freq) %>%
      rename(distance_clp = distance)
    gmp_dist <- calc_distances(gmp_freq, dest_freq) %>%
      rename(distance_gmp = distance)

    obs_df <- left_join(clp_dist, gmp_dist, by = "metric") %>%
      mutate(delta = distance_gmp - distance_clp)

    # Null hypothesis: the same pooled source mass could have been assigned to
    # CLP or GMP with probability equal to the observed global CLP fraction.
    combined_mass <- clp_freq + gmp_freq
    p_clp_global <- sum(clp_freq) / max(sum(clp_freq + gmp_freq), pseudo_count)

    null_deltas <- matrix(NA_real_, nrow = n_perms, ncol = length(metric_functions))
    colnames(null_deltas) <- names(metric_functions)

    for (perm in seq_len(n_perms)) {
      perm_clp <- rbinom(length(combined_mass),
                         size = round(combined_mass * 1e6),
                         prob = p_clp_global) / 1e6
      perm_gmp <- combined_mass - perm_clp

      perm_clp_freq <- normalize(perm_clp)
      perm_gmp_freq <- normalize(perm_gmp)

      for (metric_name in names(metric_functions)) {
        dist_fun <- metric_functions[[metric_name]]
        null_deltas[perm, metric_name] <-
          dist_fun(perm_gmp_freq, dest_freq) - dist_fun(perm_clp_freq, dest_freq)
      }
    }

    for (metric_name in names(metric_functions)) {
      obs_delta <- obs_df %>%
        filter(metric == metric_name) %>%
        pull(delta)

      p_two_sided <- mean(abs(null_deltas[, metric_name]) >= abs(obs_delta),
                          na.rm = TRUE)
      p_clp_closer <- mean(null_deltas[, metric_name] >= obs_delta,
                           na.rm = TRUE)
      p_gmp_closer <- mean(null_deltas[, metric_name] <= obs_delta,
                           na.rm = TRUE)

      metric_row <- obs_df %>%
        filter(metric == metric_name)

      cat("    ", metric_name,
          " delta:", round(obs_delta, 4),
          " | two-sided p:", round(p_two_sided, 4), "\n", sep = "")

      all_results[[length(all_results) + 1]] <- tibble(
        destination      = dest_pat,
        destination_label = dest_lab,
        bone             = bone_name,
        metric           = metric_name,
        n_clp_barcodes_used = n_clp_used,
        n_gmp_barcodes_used = n_gmp_used,
        distance_clp     = metric_row$distance_clp,
        distance_gmp     = metric_row$distance_gmp,
        delta            = metric_row$delta,
        p_two_sided      = p_two_sided,
        p_clp_closer     = p_clp_closer,
        p_gmp_closer     = p_gmp_closer
      )
    }
  }
}

results_df <- bind_rows(all_results) %>%
  group_by(destination, metric) %>%
  mutate(
    p_bh_two_sided  = p.adjust(p_two_sided, method = "BH"),
    p_bh_clp_closer = p.adjust(p_clp_closer, method = "BH"),
    p_bh_gmp_closer = p.adjust(p_gmp_closer, method = "BH")
  ) %>%
  ungroup()

write_csv(results_df,
          file.path(output_dir, "pooled_perm_CLP_vs_GMP_multi_distance_results.csv"))
cat("\nResults saved.\n")
print(results_df)

# Also save a wide table for easier review.
results_wide <- results_df %>%
  select(destination, destination_label, bone, metric,
         distance_clp, distance_gmp, delta,
         p_two_sided, p_bh_two_sided,
         p_clp_closer, p_bh_clp_closer,
         p_gmp_closer, p_bh_gmp_closer,
         n_clp_barcodes_used, n_gmp_barcodes_used) %>%
  pivot_wider(
    names_from = metric,
    values_from = c(distance_clp, distance_gmp, delta,
                    p_two_sided, p_bh_two_sided,
                    p_clp_closer, p_bh_clp_closer,
                    p_gmp_closer, p_bh_gmp_closer)
  )

write_csv(results_wide,
          file.path(output_dir, "pooled_perm_CLP_vs_GMP_multi_distance_results_wide.csv"))

# -----------------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------------

p_to_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

metric_labels <- c(
  bray_curtis    = "Bray-Curtis",
  cosine         = "Cosine",
  jensen_shannon = "Jensen-Shannon"
)

plot_df <- results_df %>%
  mutate(
    stars = p_to_stars(p_bh_two_sided),
    closer_to = if_else(delta >= 0, "CLP", "GMP"),
    metric_label = recode(metric, !!!metric_labels),
    destination_label = factor(destination_label,
                               levels = c("Blood T cells", "DN1 Thymocytes"))
  )

p_delta <- ggplot(plot_df,
                  aes(x = bone, y = delta, fill = closer_to)) +
  geom_col() +
  geom_text(aes(label = stars,
                y = delta + sign(delta) * 0.003),
            size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  facet_grid(metric_label ~ destination_label, scales = "free_x") +
  scale_fill_manual(
    values = c("CLP" = "#2196A8", "GMP" = "#E57438"),
    name = "Closer to"
  ) +
  theme_bw() +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.y     = element_text(size = 11),
    legend.position = "bottom"
  ) +
  labs(
    title    = "Pooled permutation test -- CLP versus GMP contribution to T lineage",
    subtitle = paste0("Destination vectors normalized | Source bone >= ", bone_min_umi,
                      " UMIs | ", n_perms, " permutations\n",
                      "Positive = closer to CLP | Negative = closer to GMP | Stars = BH-adjusted two-sided p"),
    x        = "Bone",
    y        = "Delta distance: distance(GMP, destination) - distance(CLP, destination)"
  )

ggsave(file.path(output_dir, "tornado_CLP_vs_GMP_multi_distance.png"),
       p_delta, width = 13, height = 9)
ggsave(file.path(output_dir, "tornado_CLP_vs_GMP_multi_distance.pdf"),
       p_delta, width = 13, height = 9)
cat("Saved: tornado_CLP_vs_GMP_multi_distance\n")

raw_plot_df <- results_df %>%
  select(destination_label, bone, metric, distance_clp, distance_gmp) %>%
  pivot_longer(cols = c(distance_clp, distance_gmp),
               names_to = "source",
               values_to = "distance") %>%
  mutate(
    source = recode(source,
                    "distance_clp" = "CLP",
                    "distance_gmp" = "GMP"),
    metric_label = recode(metric, !!!metric_labels),
    destination_label = factor(destination_label,
                               levels = c("Blood T cells", "DN1 Thymocytes"))
  )

p_raw <- ggplot(raw_plot_df,
                aes(x = bone, y = distance, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  facet_grid(metric_label ~ destination_label, scales = "free_x") +
  scale_fill_manual(
    values = c("CLP" = "#2196A8", "GMP" = "#E57438"),
    name = "Source"
  ) +
  theme_bw() +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.y     = element_text(size = 11),
    legend.position = "bottom"
  ) +
  labs(
    title = "Raw distances from CLP or GMP to each destination population",
    x = "Bone",
    y = "Distance"
  )

ggsave(file.path(output_dir, "raw_distance_CLP_vs_GMP_multi_distance.png"),
       p_raw, width = 13, height = 9)
ggsave(file.path(output_dir, "raw_distance_CLP_vs_GMP_multi_distance.pdf"),
       p_raw, width = 13, height = 9)
cat("Saved: raw_distance_CLP_vs_GMP_multi_distance\n")

cat("\nAll done! Results saved in:", output_dir, "\n")
