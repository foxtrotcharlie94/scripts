# ============================================================
# BMDM — volcano plots + up/down DE-gene barplot, from the DEG files.
#   1) volcano plots per comparison, at three significance criteria:
#        unadjusted P < 0.05, FDR < 0.25, FDR < 0.05
#      (individual files + a combined faceted figure for each criterion)
#   2) one barplot of total up / down genes per comparison
#      (filter: unadjusted PValue < 0.05)
#
# Reads <comp>/DEG_<comp>_bulk.txt (cols: gene_name, logFC, logCPM, F, PValue, FDR).
# Sign of logFC:  + = up in HB (HBvsLB)  /  + = up in lenti-pos (posvsneg)
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  have_repel <- requireNamespace("ggrepel", quietly = TRUE)
})
theme_set(theme_bw(base_size = 12))

# ---- settings ----
bmdm_dir <- "C:/Users/fc809/Downloads/BMDM"
out_dir  <- file.path(bmdm_dir, "volcano_DEcounts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COMPS  <- c("HB.neg_vs_LB.neg","HB.pos_vs_LB.pos","HB.pos_vs_HB.neg","LB.pos_vs_LB.neg")
SIGN   <- c(HB.neg_vs_LB.neg="+ = up in HB",  HB.pos_vs_LB.pos="+ = up in HB",
            HB.pos_vs_HB.neg="+ = up in pos", LB.pos_vs_LB.neg="+ = up in pos")

P_THR    <- 0.05     # unadjusted p cutoff (used for BOTH coloring and the bar counts)
LFC_THR  <- 0        # extra |logFC| cutoff for the volcano color (0 = none)
N_LABEL  <- 15       # genes to label per volcano (largest -log10P * |logFC|)

read_deg <- function(f) {
  d <- fread(file.path(bmdm_dir, f, sprintf("DEG_%s_bulk.txt", f)))
  setnames(d, 1, "ensembl")
  d <- d[, .(gene = gene_name, logFC, PValue, FDR)]
  d[, comparison := f]
  d
}
DE <- rbindlist(lapply(COMPS, read_deg))
DE[, comparison := factor(comparison, levels = COMPS)]

# ════════════════════════════════════════════════════════════════════════════
# 1) VOLCANO PLOTS
# ════════════════════════════════════════════════════════════════════════════
DE[, neglog10P := -log10(pmax(PValue, .Machine$double.xmin))]
cols_dir <- c(Up = "#C0392B", Down = "#2C6FBF", ns = "grey80")

# significance criteria -> one volcano set each
CRIT <- list(
  list(col = "PValue", cut = 0.05, name = "P0.05",   lab = "unadjusted P < 0.05"),
  list(col = "FDR",    cut = 0.25, name = "FDR0.25", lab = "FDR < 0.25"),
  list(col = "FDR",    cut = 0.05, name = "FDR0.05", lab = "FDR < 0.05")
)

volcano_one <- function(f, crit) {
  d <- copy(DE[comparison == f])
  sig <- d[[crit$col]] < crit$cut
  d[, dir := fifelse(sig & logFC >  LFC_THR, "Up",
             fifelse(sig & logFC < -LFC_THR, "Down", "ns"))]
  d[, dir := factor(dir, levels = c("Down","ns","Up"))]
  nU <- sum(d$dir == "Up"); nD <- sum(d$dir == "Down")
  # dashed line = largest P among significant genes (maps the cutoff onto the P axis)
  yb <- if (any(sig)) -log10(max(d$PValue[sig], na.rm = TRUE)) else NA_real_
  lab <- d[sig][order(-(neglog10P * abs(logFC)))][seq_len(min(N_LABEL, .N))]
  p <- ggplot(d, aes(logFC, neglog10P, color = dir)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = cols_dir, name = NULL, drop = FALSE) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    labs(title = f,
         subtitle = sprintf("%s | %s | up %d, down %d", SIGN[[f]], crit$lab, nU, nD),
         x = "log2 fold-change", y = "-log10(P-value)")
  if (!is.na(yb)) p <- p + geom_hline(yintercept = yb, linetype = "dashed", color = "grey50")
  if (have_repel && nrow(lab) > 0)
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene),
              size = 3, color = "grey20", max.overlaps = 40, segment.color = "grey70")
  p
}

for (crit in CRIT) {
  vlist <- lapply(COMPS, volcano_one, crit = crit)
  for (i in seq_along(COMPS))                     # individual volcano files
    ggsave(file.path(out_dir, sprintf("volcano_%s_%s.png", COMPS[i], crit$name)),
           vlist[[i]], width = 7, height = 6, dpi = 150)
  combined <- wrap_plots(vlist, ncol = 2) + plot_layout(guides = "collect") +
    plot_annotation(title = sprintf("BMDM volcano plots — colored by %s", crit$lab))
  ggsave(file.path(out_dir, sprintf("volcano_all_%s.png", crit$name)),
         combined, width = 13, height = 11, dpi = 150)
  cat(sprintf("volcano set '%s' written\n", crit$lab))
}

# ════════════════════════════════════════════════════════════════════════════
# 2) UP / DOWN GENE COUNTS (unadjusted P < 0.05)
# ════════════════════════════════════════════════════════════════════════════
counts <- DE[PValue < P_THR, .(
  Up   = sum(logFC > 0),
  Down = sum(logFC < 0)), by = comparison]
fwrite(counts, file.path(out_dir, "DE_counts_P0.05.csv"))
cat("Up/Down genes at P <", P_THR, ":\n"); print(counts)

cnt_long <- melt(counts, id.vars = "comparison", variable.name = "direction", value.name = "n")
cnt_long[, direction := factor(direction, levels = c("Up","Down"))]
cnt_long[,