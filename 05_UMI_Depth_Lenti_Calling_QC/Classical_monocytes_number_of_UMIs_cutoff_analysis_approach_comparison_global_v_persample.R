## ============================================================
##  Lenti detection threshold comparison — LT-HSCs
##  Approach A: top 50% globally (pool all cells, rank globally)
##  Approach B: top 50% per sample (rank within each sample)
##  Output: detection curves + printed thresholds for both
## ============================================================

library(data.table)
library(tidyverse)
library(MASS)
library(ggplot2)

setwd("C:/Users/fc809/Downloads/")

MATRIX_FILE <- "LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt"
DETECT_PROB <- 0.80
TOP_FRAC    <- 0.50

## ── 1. Read matrix ────────────────────────────────────────────────────────────
message("Reading matrix...")
mat_raw    <- as.data.frame(fread(MATRIX_FILE))
gene_names <- mat_raw[[1]]
mat        <- as.matrix(mat_raw[, -1])
rownames(mat) <- gene_names
barcodes   <- colnames(mat)
message("  ", nrow(mat), " genes x ", ncol(mat), " cells")

lenti_row  <- which(grepl("^[Ll]enti$", gene_names))[1]
if (is.na(lenti_row)) lenti_row <- which(grepl("[Ll]enti", gene_names))[1]
message("  Lenti row: '", gene_names[lenti_row], "'")

numi_vec   <- colSums(mat)
lenti_umis <- mat[lenti_row, ]
lenti_freq <- lenti_umis / numi_vec

df <- data.frame(
  barcode    = barcodes,
  sample     = sub("_.*", "", barcodes),
  nUMI       = numi_vec,
  lenti_UMIs = as.integer(lenti_umis),
  lenti_freq = lenti_freq,
  stringsAsFactors = FALSE
)

numi_range <- seq(1000, 150000, by = 100)

## ── Helper: fit NB and return detection curve + threshold ─────────────────────
fit_nb_curve <- function(training_df, label) {
  training_df$log_nUMI <- log(training_df$nUMI)
  fit <- tryCatch(
    glm.nb(lenti_UMIs ~ log_nUMI, data = training_df, link = log),
    error = function(e) {
      message("  glm.nb failed for '", label, "', falling back to Poisson")
      glm(lenti_UMIs ~ log_nUMI, data = training_df, family = poisson(link = "log"))
    }
  )
  b0    <- coef(fit)[["(Intercept)"]]
  b1    <- coef(fit)[["log_nUMI"]]
  theta <- if (inherits(fit, "negbin")) fit$theta else Inf
  
  probs <- sapply(numi_range, function(n) {
    mu <- exp(b0 + b1 * log(n))
    if (is.finite(theta)) 1 - (theta / (theta + mu))^theta
    else 1 - exp(-mu)
  })
  
  thresh <- approx(probs, numi_range, xout = DETECT_PROB)$y
  list(probs = probs, threshold = thresh, b0 = b0, b1 = b1, theta = theta)
}

## ── 2. Approach A: global top 50% ─────────────────────────────────────────────
message("\n=== Approach A: global top 50% ===")
n_pos_A <- ceiling(nrow(df) * TOP_FRAC)
df_A    <- df[order(df$lenti_freq, decreasing = TRUE), ][seq_len(n_pos_A), ]
message("  Training cells: ", nrow(df_A))
res_A   <- fit_nb_curve(df_A, "global")
message("  Threshold: ", round(res_A$threshold), " UMIs")

## ── 3. Approach B: per-sample top 50%, then pool ─────────────────────────────
message("\n=== Approach B: per-sample top 50% ===")
samples <- sort(unique(df$sample))
df_B_list <- lapply(samples, function(samp) {
  ds      <- df[df$sample == samp, ]
  n_pos_s <- ceiling(nrow(ds) * TOP_FRAC)
  ds_top  <- ds[order(ds$lenti_freq, decreasing = TRUE), ][seq_len(n_pos_s), ]
  message("  ", samp, ": ", nrow(ds), " cells → ", nrow(ds_top), " training")
  ds_top
})
df_B <- bind_rows(df_B_list)
message("  Total training cells: ", nrow(df_B))
res_B <- fit_nb_curve(df_B, "per-sample pooled")
message("  Threshold: ", round(res_B$threshold), " UMIs")

## ── 4. Summary ────────────────────────────────────────────────────────────────
message("\n=== COMPARISON ===")
message(sprintf("  %-35s  %8s", "Approach", "Threshold"))
message(sprintf("  %-35s  %8d", "A: global top 50%", round(res_A$threshold)))
message(sprintf("  %-35s  %8d", "B: per-sample top 50% (pooled)", round(res_B$threshold)))
message(sprintf("  %-35s  %8s", "Difference (B - A)",
                paste0(ifelse(res_B$threshold - res_A$threshold > 0, "+", ""),
                       round(res_B$threshold - res_A$threshold))))

## ── 5. Plot: two curves only ──────────────────────────────────────────────────
curve_df <- bind_rows(
  data.frame(nUMI = numi_range, prob = res_A$probs,
             approach = paste0("A: global top 50%\n(threshold = ",
                               round(res_A$threshold), ")")),
  data.frame(nUMI = numi_range, prob = res_B$probs,
             approach = paste0("B: per-sample top 50%\n(threshold = ",
                               round(res_B$threshold), ")"))
)

cols <- c("A: global top 50%\n(threshold = " = "#2C3E50",
          "B: per-sample top 50%\n(threshold = " = "#C0392B")

# Use a 2-colour palette keyed on approach
approach_levels <- unique(curve_df$approach)
approach_cols   <- setNames(c("#2C3E50", "#C0392B"), approach_levels)

p <- ggplot(curve_df, aes(x = nUMI, y = prob, color = approach)) +
  geom_line(linewidth = 1.4) +
  geom_hline(yintercept = DETECT_PROB, linetype = "dashed",
             color = "grey40", linewidth = 0.8) +
  geom_vline(xintercept = res_A$threshold, linetype = "dotted",
             color = "#2C3E50", linewidth = 0.8) +
  geom_vline(xintercept = res_B$threshold, linetype = "dotted",
             color = "#C0392B", linewidth = 0.8) +
  annotate("text", x = res_A$threshold, y = 0.04,
           label = round(res_A$threshold),
           hjust = 0.5, size = 4, color = "#2C3E50", fontface = "bold") +
  annotate("text", x = res_B$threshold, y = 0.10,
           label = round(res_B$threshold),
           hjust = 0.5, size = 4, color = "#C0392B", fontface = "bold") +
  scale_color_manual(values = approach_cols, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title    = "Lenti Detection Threshold Comparison — LT-HSCs",
    subtitle = paste0("Both curves fitted on top-50% LentiPos cells (NB model)\n",
                      "Dashed grey = 80% detection threshold"),
    x = "Total UMI Count per Cell",
    y = "P(detect \u2265 1 Lenti UMI)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 16),
    plot.subtitle   = element_text(size = 12, color = "grey40"),
    legend.position = "bottom",
    legend.text     = element_text(size = 12)
  )

ggsave("Lenti_threshold_comparison_LTHSCs.pdf", p, width = 9, height = 6)
ggsave("Lenti_threshold_comparison_LTHSCs.png", p, width = 9, height = 6, dpi = 180)

message("\nPlot saved.")
