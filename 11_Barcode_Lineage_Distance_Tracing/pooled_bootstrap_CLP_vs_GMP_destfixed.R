# =============================================================================
# Pooled bootstrap — CLP vs GMP contribution to T lineage (destination FIXED)
# Metrics: Bray-Curtis, Cosine, Jensen-Shannon
#
# Question: holding the observed destination (Blood T / Thymus DN1) fixed, is its
# composition reliably closer to CLP than to GMP, given how well we sampled the
# two SOURCES?  We challenge only the source-side sampling, not the destination
# (the destination is the fixed target we are attributing flux INTO).
#
# Method: treat each pooled source as a multinomial at its OBSERVED read depth
# (CLP at N_CLP UMIs, GMP at N_GMP UMIs), draw bootstrap replicates, renormalize,
# and recompute delta with the destination held fixed.  Report the observed delta
# and a percentile CI; the ranking is reliable when the CI excludes 0.
#
#   delta = distance(GMP, destination) - distance(CLP, destination)
#   delta > 0  =>  destination closer to CLP
#   delta < 0  =>  destination closer to GMP
#
# What this CI captures: read-sampling noise in the two sources at the pooled
# depth.  It does NOT capture biological variation across mice — that is the job
# of the per-mouse companion analysis (mice as replicates).  Uses UNFILTERED
# barcodes, source-side per-sample threshold >= bone_min_umi, pooled across mice.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# SETTINGS  (paths as in the original script — edit for your machine)
# -----------------------------------------------------------------------------
base_dir   <- "/Users/karingustafsson/Documents/Jobb/2026/Results/Sternum/DARLIN"
output_dir <- file.path(base_dir, "pooled_bootstrap_CLP_vs_GMP_destfixed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

bones_of_interest  <- c("Femur", "Sternum", "Pelvis", "Spine")
source_populations <- c("CLP", "GMP")
bone_min_umi       <- 2
n_boot             <- 2000     # bootstrap replicates (2000 is ample for a 95% CI)
seed               <- 1
ci_level           <- 0.95

destinations <- tribble(
  ~dest_pattern,     ~dest_label,
  "Blood_T_cells",   "Blood T cells",
  "Thymus_DN1",      "DN1 Thymocytes"
)

# -----------------------------------------------------------------------------
# LOAD UNFILTERED BARCODE COUNTS
# (C and T concatenated and normalized together => C/T contribute in proportion
#  to their raw UMI depth.  To force a balanced 50/50 C:T instead, normalize each
#  region separately before combining — left as the original depth-weighted form.)
# -----------------------------------------------------------------------------
cat("Reading unfiltered C and T tables...\n")
c_wide_raw <- read_tsv(file.path(base_dir, "combined_barcodes_C_regions.tsv"), show_col_types = FALSE)
t_wide_raw <- read_tsv(file.path(base_dir, "combined_barcodes_T_regions.tsv"), show_col_types = FALSE)

strip_suffix <- function(x) gsub("_(CA|TA)$", "", gsub("_S\\d+$", "", x))
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
  for (col in setdiff(all_cols, colnames(df)[-1])) df[[col]] <- 0
  df %>% select(barcode, all_of(all_cols))
}

raw_counts <- bind_rows(add_missing(c_counts, all_samples),
                        add_missing(t_counts, all_samples))
cat("Wide table:", nrow(raw_counts), "barcodes x", ncol(raw_counts), "columns\n")

# -----------------------------------------------------------------------------
# SAMPLE METADATA
# -----------------------------------------------------------------------------
sample_meta <- tibble(sample = all_samples) %>%
  mutate(
    mouse      = map_chr(sample, ~ str_split(.x, "_")[[1]][1]),
    population = map_chr(sample, ~ str_split(.x, "_")[[1]][2]),
    bone       = map_chr(sample, ~ str_split(.x, "_")[[1]][3])
  ) %>%
  filter(!is.na(suppressWarnings(as.integer(mouse)))) %>%
  mutate(bone = str_replace(bone, "Pelivis", "Pelvis"))
cat("Mice found:", paste(unique(sample_meta$mouse), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# DISTANCE FUNCTIONS  (matrix form: F = barcodes x replicates, dest = fixed vector)
# -----------------------------------------------------------------------------
bc_mat <- function(F, dest) colSums(abs(F - dest)) / colSums(F + dest)           # dest recycled down columns
cos_mat <- function(F, dest) as.numeric(1 - crossprod(dest, F) /
                                          (sqrt(colSums(F^2)) * sqrt(sum(dest^2))))
js_one <- function(p, q) {
  if (sum(p) == 0 || sum(q) == 0) return(NA_real_)
  p <- p / sum(p); q <- q / sum(q); m <- (p + q) / 2
  kl <- function(x, y) { i <- x > 0 & y > 0; if (!any(i)) 0 else sum(x[i] * log2(x[i] / y[i])) }
  sqrt(max(0, 0.5 * kl(p, m) + 0.5 * kl(q, m)))
}
js_mat <- function(F, dest) apply(F, 2, js_one, q = dest)

dist_mat_funs <- list(bray_curtis = bc_mat, cosine = cos_mat, jensen_shannon = js_mat)

# -----------------------------------------------------------------------------
# POOL UMIs ACROSS MICE  (per-sample threshold applied before pooling)
# -----------------------------------------------------------------------------
pool_source <- function(population_name, bone_name) {
  ss <- sample_meta %>% filter(population == population_name, bone == bone_name) %>% pull(sample)
  if (length(ss) == 0) return(rep(0, nrow(raw_counts)))
  mat <- as.matrix(select(raw_counts, all_of(ss)))
  mat[mat < bone_min_umi] <- 0
  rowSums(mat)
}
pool_dest <- function(dest_pattern) {
  ds <- sample_meta %>% filter(str_detect(sample, dest_pattern)) %>% pull(sample)
  if (length(ds) == 0) return(rep(0, nrow(raw_counts)))
  rowSums(as.matrix(select(raw_counts, all_of(ds))))
}

# -----------------------------------------------------------------------------
# MAIN: observed delta + source-side bootstrap CI, destination held fixed
# -----------------------------------------------------------------------------
set.seed(seed)
alpha <- (1 - ci_level) / 2
all_results <- list()

for (di in seq_len(nrow(destinations))) {
  dest_pat <- destinations$dest_pattern[di]
  dest_lab <- destinations$dest_label[di]
  cat("\n=== Destination:", dest_lab, "===\n")

  dest_full <- pool_dest(dest_pat)
  if (sum(dest_full) == 0) { cat("  No destination UMIs -- skipping\n"); next }
  cat("  Destination barcodes:", sum(dest_full > 0), "| total UMIs:", sum(dest_full), "\n")

  for (bone_name in bones_of_interest) {
    clp_full <- pool_source("CLP", bone_name)
    gmp_full <- pool_source("GMP", bone_name)
    n_clp_used <- sum(clp_full > 0); n_gmp_used <- sum(gmp_full > 0)
    if (n_clp_used == 0 || n_gmp_used == 0) {
      cat("  Bone:", bone_name, "-- missing CLP or GMP barcodes >=", bone_min_umi, "UMIs -- skipping\n"); next
    }

    # Restrict to barcodes nonzero in CLP, GMP, or destination (exact: all-zero
    # barcodes contribute nothing to any of the three distances).
    active   <- which(clp_full > 0 | gmp_full > 0 | dest_full > 0)
    clp_c    <- clp_full[active]; gmp_c <- gmp_full[active]; dst_c <- dest_full[active]
    N_clp    <- round(sum(clp_c)); N_gmp <- round(sum(gmp_c)); N_dst <- sum(dst_c)
    clp_p    <- clp_c / sum(clp_c); gmp_p <- gmp_c / sum(gmp_c)
    dest_freq <- dst_c / N_dst                     # FIXED across all replicates

    cat("  Bone:", bone_name,
        "| N_CLP:", N_clp, "| N_GMP:", N_gmp, "| N_dest:", N_dst,
        "| active barcodes:", length(active), "\n")

    # ---- bootstrap source replicates at real depth ----
    clp_boot <- rmultinom(n_boot, N_clp, clp_p)
    gmp_boot <- rmultinom(n_boot, N_gmp, gmp_p)
    clp_bf <- sweep(clp_boot, 2, colSums(clp_boot), "/")
    gmp_bf <- sweep(gmp_boot, 2, colSums(gmp_boot), "/")

    for (metric_name in names(dist_mat_funs)) {
      fn <- dist_mat_funs[[metric_name]]

      # observed (plug-in)
      d_clp_obs <- fn(matrix(clp_p, ncol = 1), dest_freq)
      d_gmp_obs <- fn(matrix(gmp_p, ncol = 1), dest_freq)
      delta_obs <- d_gmp_obs - d_clp_obs

      # bootstrap distributions
      dC <- fn(clp_bf, dest_freq)
      dG <- fn(gmp_bf, dest_freq)
      delta_b <- dG - dC

      ci_d  <- quantile(delta_b, c(alpha, 1 - alpha), na.rm = TRUE)
      ciC   <- quantile(dC,      c(alpha, 1 - alpha), na.rm = TRUE)
      ciG   <- quantile(dG,      c(alpha, 1 - alpha), na.rm = TRUE)
      p_two <- min(1, 2 * min(mean(delta_b <= 0, na.rm = TRUE),
                              mean(delta_b >= 0, na.rm = TRUE)))
      prob_clp_closer <- mean(delta_b > 0, na.rm = TRUE)
      ci_excludes_0   <- (ci_d[1] > 0) || (ci_d[2] < 0)

      cat("    ", metric_name,
          " delta:", round(delta_obs, 4),
          " CI[", round(ci_d[1], 4), ",", round(ci_d[2], 4), "]",
          if (ci_excludes_0) " *" else "", "\n", sep = "")

      all_results[[length(all_results) + 1]] <- tibble(
        destination = dest_pat, destination_label = dest_lab, bone = bone_name,
        metric = metric_name,
        N_clp = N_clp, N_gmp = N_gmp, N_dest = N_dst,
        n_clp_barcodes_used = n_clp_used, n_gmp_barcodes_used = n_gmp_used,
        distance_clp = d_clp_obs, distance_clp_lo = ciC[1], distance_clp_hi = ciC[2],
        distance_gmp = d_gmp_obs, distance_gmp_lo = ciG[1], distance_gmp_hi = ciG[2],
        delta = delta_obs, delta_ci_lo = ci_d[1], delta_ci_hi = ci_d[2],
        boot_p_two_sided = p_two, prob_clp_closer = prob_clp_closer,
        ci_excludes_0 = ci_excludes_0
      )
    }
  }
}

results_df <- bind_rows(all_results) %>%
  group_by(destination, metric) %>%
  mutate(boot_p_bh = p.adjust(boot_p_two_sided, method = "BH")) %>%
  ungroup()

write_csv(results_df, file.path(output_dir, "pooled_bootstrap_CLP_vs_GMP_destfixed_results.csv"))
cat("\nResults saved.\n"); print(results_df)

# -----------------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------------
metric_labels <- c(bray_curtis = "Bray-Curtis", cosine = "Cosine", jensen_shannon = "Jensen-Shannon")

plot_df <- results_df %>%
  mutate(
    closer_to = if_else(delta >= 0, "CLP", "GMP"),
    metric_label = recode(metric, !!!metric_labels),
    destination_label = factor(destination_label, levels = c("Blood T cells", "DN1 Thymocytes")),
    sig = if_else(ci_excludes_0, "*", "")
  )

p_delta <- ggplot(plot_df, aes(x = bone, y = delta, fill = closer_to)) +
  geom_col() +
  geom_errorbar(aes(ymin = delta_ci_lo, ymax = delta_ci_hi), width = 0.25, color = "gray20") +
  geom_text(aes(label = sig, y = delta_ci_hi + sign(delta) * 0.004), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  facet_grid(metric_label ~ destination_label, scales = "free_x") +
  scale_fill_manual(values = c(CLP = "#2196A8", GMP = "#E57438"), name = "Closer to") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 11), legend.position = "bottom") +
  labs(title = "CLP vs GMP contribution to T lineage — destination fixed, source bootstrap",
       subtitle = paste0("Source bone >= ", bone_min_umi, " UMIs | ", n_boot,
                         " bootstraps | ", round(ci_level * 100), "% percentile CI\n",
                         "Positive = closer to CLP | * = CI excludes 0"),
       x = "Bone",
       y = "Delta: distance(GMP, dest) - distance(CLP, dest)")

ggsave(file.path(output_dir, "tornado_CLP_vs_GMP_destfixed.png"), p_delta, width = 13, height = 9)
ggsave(file.path(output_dir, "tornado_CLP_vs_GMP_destfixed.pdf"), p_delta, width = 13, height = 9)

raw_plot_df <- bind_rows(
  results_df %>% transmute(destination_label, bone, metric, source = "CLP",
                           distance = distance_clp, lo = distance_clp_lo, hi = distance_clp_hi),
  results_df %>% transmute(destination_label, bone, metric, source = "GMP",
                           distance = distance_gmp, lo = distance_gmp_lo, hi = distance_gmp_hi)
) %>%
  mutate(metric_label = recode(metric, !!!metric_labels),
         destination_label = factor(destination_label, levels = c("Blood T cells", "DN1 Thymocytes")))

p_raw <- ggplot(raw_plot_df, aes(x = bone, y = distance, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8),
                width = 0.25, color = "gray20") +
  coord_flip() +
  facet_grid(metric_label ~ destination_label, scales = "free_x") +
  scale_fill_manual(values = c(CLP = "#2196A8", GMP = "#E57438"), name = "Source") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 11), legend.position = "bottom") +
  labs(title = "Raw distances from CLP / GMP to each destination (with source-bootstrap CIs)",
       x = "Bone", y = "Distance")

ggsave(file.path(output_dir, "raw_distance_CLP_vs_GMP_destfixed.png"), p_raw, width = 13, height = 9)
ggsave(file.path(output_dir, "raw_distance_CLP_vs_GMP_destfixed.pdf"), p_raw, width = 13, height = 9)

cat("\nAll done! Results saved in:", output_dir, "\n")
