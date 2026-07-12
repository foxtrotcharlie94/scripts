# ============================================================
# Shared-axis pathway figure: HB-vs-LB (burden) is the dominant axis of
# pathway change, and LentiPos-vs-LentiNeg (lenti) is a partial move along
# the SAME axis (same direction, smaller magnitude).
#
# Works entirely from the GSEA_*.csv files in pop_dir (fast, CSV-only).
#
#   Panel A  burden NES (x) vs lenti NES (y), one point per pathway.
#            TLS-through-origin slope = "dose" (lenti as a fraction of a
#            full burden step), with a pathway-bootstrap CI. Quadrants +
#            labels show the shared mechanism.
#   Panel B  per-contrast decomposition: bar height = NES-vector length
#            (size of the pathway response = dominance); dark portion =
#            projection onto the burden axis (how much is the burden program).
#
# Headline = Hallmark, averaged axes, logFC. GO_BP / logFCxP / per-condition
# versions are written as supplements.
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
have_repel <- requireNamespace("ggrepel", quietly = TRUE)

# ── Settings ──────────────────────────────────────────────────────────────────
pop_dir      <- "C:/Users/fc809/Downloads/pseudobulk_DE_GSEA/ClassicalMonocytes"

burden_comps <- c("LentiNeg_HBvsLB", "LentiPos_HBvsLB")   # +NES = up in HB
lenti_comps  <- c("HB_PosVsNeg",     "LB_PosVsNeg")        # +NES = up in LentiPos
all_comps    <- c(burden_comps, lenti_comps)

collections <- c("Hallmark", "GO_BP")    # Hallmark = headline, GO_BP = supplement
rankings    <- c("logFC", "logFCxP")     # logFC = headline
head_coll <- "Hallmark"; head_rank <- "logFC"

sig_strong <- 0.05    # INCLUSION part 1: FDR < this in >=1 contrast
sig_loose  <- 0.10    # INCLUSION part 2: FDR < this in ALL contrasts (consistency)
n_label  <- 10        # pathways to label in the scatter (largest |burden|)
n_boot   <- 2000

out_dir <- file.path(pop_dir, "_shared_axis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────
gsea_path <- function(comp, coll, rk)
  file.path(pop_dir, comp, sprintf("GSEA_%s_%s_rankBy_%s.csv", coll, comp, rk))

load_one <- function(comp, coll, rk) {
  fp <- gsea_path(comp, coll, rk)
  if (!file.exists(fp)) { warning("Missing: ", fp); return(NULL) }
  d <- fread(fp)[, .(pathway, NES, pval, padj)]
  setnames(d, c("NES","pval","padj"),
           c(paste0("NES_", comp), paste0("pval_", comp), paste0("padj_", comp))); d
}

clean_name <- function(p) gsub("_", " ", sub("^GOBP_", "", sub("^HALLMARK_", "", p)))

tls0 <- function(x, y) {            # total-least-squares slope through origin
  ev <- eigen(crossprod(cbind(x, y)))$vectors[, 1]; ev[2] / ev[1]
}

# ── Accumulators ──────────────────────────────────────────────────────────────
dose_rows <- list(); decomp_rows <- list()

for (coll in collections) for (rk in rankings) {
  parts <- lapply(all_comps, load_one, coll = coll, rk = rk)
  if (any(sapply(parts, is.null))) { cat(sprintf("[%s|%s] missing files, skipping\n", coll, rk)); next }
  m <- Reduce(function(a, b) merge(a, b, by = "pathway"), parts)
  m <- m[is.finite(rowSums(as.matrix(m[, grep("^NES_", names(m)), with = FALSE])))]
  pcols <- paste0("padj_", all_comps)
  P <- as.matrix(m[, ..pcols])
  m <- m[(rowSums(P < sig_strong, na.rm = TRUE) >= 1) &                       # FDR<0.05 in >=1
         (rowSums(P < sig_loose,  na.rm = TRUE) == length(all_comps))]        # AND FDR<0.10 in all
  if (nrow(m) < 8) { cat(sprintf("[%s|%s] too few pathways after FDR<%.2g(any)+<%.2g(all) filter\n", coll, rk, sig_strong, sig_loose)); next }

  NES <- as.matrix(m[, paste0("NES_", all_comps), with = FALSE]); rownames(NES) <- m$pathway
  burden <- rowMeans(NES[, paste0("NES_", burden_comps)])
  lenti  <- rowMeans(NES[, paste0("NES_", lenti_comps)])

  # dose: lenti ~ slope * burden (through origin, TLS)
  slope <- tls0(burden, lenti)
  set.seed(1)
  sb <- replicate(n_boot, { i <- sample(length(burden), replace = TRUE); tls0(burden[i], lenti[i]) })
  ci <- quantile(sb, c(.025, .975), na.rm = TRUE)
  pear <- cor(burden, lenti); spear <- cor(burden, lenti, method = "spearman")
  mag_ratio <- sqrt(sum(burden^2) / sum(lenti^2))   # ||burden|| / ||lenti||
  dose_rows[[length(dose_rows)+1]] <- data.table(
    collection = coll, ranking = rk, n_pathways = nrow(m),
    dose = slope, ci_lo = ci[1], ci_hi = ci[2],
    pearson = pear, spearman = spear, magnitude_ratio = mag_ratio)

  # decomposition: project each contrast onto the unit burden axis
  u <- burden / sqrt(sum(burden^2))
  padj_m <- as.matrix(m[, paste0("padj_", all_comps), with = FALSE])
  for (j in seq_along(all_comps)) {
    v <- NES[, paste0("NES_", all_comps[j])]
    total <- sqrt(sum(v^2)); proj <- sum(v * u)
    decomp_rows[[length(decomp_rows)+1]] <- data.table(
      collection = coll, ranking = rk, contrast = all_comps[j],
      axis = if (all_comps[j] %in% burden_comps) "burden" else "lenti",
      total = total, projection = proj, cosine = proj / total,
      n_sig = sum(padj_m[, j] < sig_strong, na.rm = TRUE))
  }

  # ── Panel A : scatter ───────────────────────────────────────────────────────
  d <- data.table(pathway = rownames(NES), burden = burden, lenti = lenti,
                  sig = (padj_m[, paste0("padj_", burden_comps[1])] < sig_strong) |
                        (padj_m[, paste0("padj_", lenti_comps[1])]  < sig_strong))
  d[, quad := fifelse(burden > 0 & lenti > 0, "both up",
               fifelse(burden < 0 & lenti < 0, "both down", "discordant"))]
  lim <- max(abs(c(d$burden, d$lenti))) * 1.05
  lab <- d[order(-abs(burden))][seq_len(min(n_label, .N))]
  quad_cols <- c("both up" = "#C44E52", "both down" = "#4C72B0", "discordant" = "grey70")

  pA <- ggplot(d, aes(burden, lenti)) +
    annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim, fill = "#C44E52", alpha = 0.05) +
    annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0, fill = "#4C72B0", alpha = 0.05) +
    geom_hline(yintercept = 0, color = "grey75") + geom_vline(xintercept = 0, color = "grey75") +
    geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey60") +
    geom_abline(slope = slope, intercept = 0, color = "black", linewidth = 0.8) +
    geom_point(aes(color = quad, alpha = sig), size = 2) +
    scale_color_manual(values = quad_cols, name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.35), guide = "none") +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    annotate("text", x = -lim, y = lim,
             label = sprintf("dose = %.2f  (95%% CI %.2f-%.2f)\nPearson r = %.2f\nburden:lenti magnitude = %.1fx",
                             slope, ci[1], ci[2], pear, mag_ratio),
             hjust = 0, vjust = 1, size = 3.2) +
    labs(title = sprintf("Lenti effect is a partial move along the burden axis (%s, rankBy %s)", coll, rk),
         subtitle = sprintf("point = pathway FDR<%.2g in \u22651 & FDR<%.2g in all (n=%d) | x = HB-vs-LB NES (avg) | y = LentiPos-vs-LentiNeg NES (avg)", sig_strong, sig_loose, nrow(d)),
         x = "burden NES  (+ = up in HB)", y = "lenti NES  (+ = up in LentiPos)") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  pA <- if (have_repel)
    pA + ggrepel::geom_text_repel(data = lab, aes(label = clean_name(pathway)), size = 2.6, max.overlaps = 20)
  else
    pA + geom_text(data = lab, aes(label = clean_name(pathway)), size = 2.4, vjust = -0.7, check_overlap = TRUE)

  # ── Panel B : decomposition ─────────────────────────────────────────────────
  db <- rbindlist(decomp_rows)[collection == coll & ranking == rk]
  db[, contrast := factor(contrast, levels = all_comps)]
  bar <- rbindlist(list(
    db[, .(contrast, axis, part = "along burden axis", value = pmax(projection, 0))],
    db[, .(contrast, axis, part = "off-axis",          value = pmax(total - pmax(projection, 0), 0))]
  ))
  bar[, part := factor(part, levels = c("off-axis", "along burden axis"))]
  pB <- ggplot(bar, aes(contrast, value, fill = part)) +
    geom_col(width = 0.7) +
    geom_text(data = db, aes(contrast, total, label = sprintf("%.0f%%", 100 * cosine)),
              inherit.aes = FALSE, vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("along burden axis" = "#2C7A6F", "off-axis" = "grey80"), name = NULL) +
    labs(title = "Size of pathway response, split by burden alignment",
         subtitle = "bar = ||NES|| (response size) | dark = projection onto burden axis | label = % aligned (cosine)",
         x = NULL, y = "pathway-response magnitude  ||NES||") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)))

  tag <- sprintf("%s_rankBy_%s", coll, rk)
  combined <- (pA | pB) + plot_layout(widths = c(1.15, 1)) +
    plot_annotation(title = sprintf("Shared-axis summary — %s | %s | rankBy %s%s",
                                    basename(pop_dir), coll, rk,
                                    if (coll == head_coll && rk == head_rank) "   [HEADLINE]" else "   [supplement]"))
  ggsave(file.path(out_dir, sprintf("shared_axis_%s.png", tag)), combined,
         width = 14, height = 7, dpi = 150)

  # ── Supplement: per-condition (keep lenti contrasts separate) ───────────────
  sep <- rbindlist(lapply(lenti_comps, function(lc)
    data.table(lenti_contrast = lc, pathway = rownames(NES),
               burden = burden, lenti = NES[, paste0("NES_", lc)])))
  sl_txt <- sep[, .(slope = tls0(burden, lenti), r = cor(burden, lenti)), by = lenti_contrast]
  pS <- ggplot(sep, aes(burden, lenti)) +
    geom_hline(yintercept = 0, color = "grey75") + geom_vline(xintercept = 0, color = "grey75") +
    geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey60") +
    geom_point(color = "#4C72B0", alpha = 0.6, size = 1.8) +
    geom_abline(data = sl_txt, aes(slope = slope, intercept = 0), color = "black", linewidth = 0.7) +
    geom_text(data = sl_txt, aes(label = sprintf("dose=%.2f, r=%.2f", slope, r)),
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 3) +
    facet_wrap(~ lenti_contrast) + coord_equal() +
    labs(title = sprintf("Lenti effect vs burden axis, by condition (%s, rankBy %s)", coll, rk),
         subtitle = "consistency check: same direction in both HB and LB",
         x = "burden NES (avg, + = up in HB)", y = "lenti NES (+ = up in LentiPos)") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, sprintf("shared_axis_byCondition_%s.png", tag)), pS,
         width = 11, height = 5.5, dpi = 150)

  cat(sprintf("[%s|%s] %d pathways | dose=%.2f (%.2f-%.2f) | r=%.2f | burden:lenti=%.1fx\n",
              coll, rk, nrow(m), slope, ci[1], ci[2], pear, mag_ratio))
}

fwrite(rbindlist(dose_rows),   file.path(out_dir, "shared_axis_dose.csv"))
fwrite(rbindlist(decomp_rows), file.path(out_dir, "shared_axis_decomposition.csv"))
cat("\nDone. Headline: shared_axis_Hallmark_rankBy_logFC.png\nAll outputs in:", out_dir, "\n")
