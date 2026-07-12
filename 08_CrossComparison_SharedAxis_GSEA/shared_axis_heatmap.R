# ============================================================
# Cross-contrast pathway heatmap.
#
# Rows  = pathways significant (FDR < sig) in AT LEAST ONE of the four contrasts.
#         Ordered by a cross-contrast aggregate score = mean over contrasts of
#         NES x -log10(FDR)  (up block on top, down block below).
# Cols  = all four contrasts.
# Fill  = NES (sign = direction, intensity = magnitude).
# Stars = FDR significance (*** <.001, ** <.01, * <.05, . <.10).
#
# Reads the GSEA_*.csv files in pop_dir (fast, CSV-only).
# ============================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# ── Settings ──────────────────────────────────────────────────────────────────
pop_dir  <- "C:/Users/fc809/Downloads/pseudobulk_DE_GSEA/ClassicalMonocytes"
comps    <- c("LentiNeg_HBvsLB", "LentiPos_HBvsLB", "HB_PosVsNeg", "LB_PosVsNeg")  # column order
burden_n <- 2          # first N columns are the burden contrasts (draw a separator after)

collections <- c("Hallmark", "GO_BP")
rankings    <- c("logFC", "logFCxP")

sig          <- 0.05   # row set: FDR < this in >=1 contrast
top_n_per_dir <- 80    # cap rows per direction for readability (Inf = no cap)
grey_nes     <- 1.5    # grey out tiles that are non-sig (FDR>=sig) AND |NES| < this

out_dir <- file.path(pop_dir, "_shared_axis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────
gsea_path <- function(comp, coll, rk)
  file.path(pop_dir, comp, sprintf("GSEA_%s_%s_rankBy_%s.csv", coll, comp, rk))
load_gsea <- function(comp, coll, rk) {
  fp <- gsea_path(comp, coll, rk); if (!file.exists(fp)) return(NULL)
  fread(fp)[, .(pathway, NES, padj)]
}
clean_name <- function(p) gsub("_", " ", sub("^GOBP_", "", sub("^HALLMARK_", "", p)))
stars <- function(p) fifelse(is.na(p), "", fifelse(p < .001, "***",
                                                   fifelse(p < .01, "**", fifelse(p < .05, "*", fifelse(p < .10, "\u00b7", "")))))

for (coll in collections) for (rk in rankings) {
  tabs <- lapply(comps, load_gsea, coll = coll, rk = rk); names(tabs) <- comps
  if (any(sapply(tabs, is.null))) { cat(sprintf("[%s|%s] missing GSEA file(s), skipping\n", coll, rk)); next }

  # row set: union of pathways FDR<sig in any contrast
  sig_any <- unique(unlist(lapply(tabs, function(g) g[!is.na(padj) & padj < sig, pathway])))
  if (length(sig_any) == 0) { cat(sprintf("[%s|%s] no sig pathways anywhere\n", coll, rk)); next }

  # wide NES / padj over the union, then aggregate ordering score
  NESm  <- sapply(comps, function(cp) tabs[[cp]]$NES [match(sig_any, tabs[[cp]]$pathway)])
  padjm <- sapply(comps, function(cp) tabs[[cp]]$padj[match(sig_any, tabs[[cp]]$pathway)])
  agg <- rowMeans(NESm * -log10(pmax(padjm, 1e-300)), na.rm = TRUE)   # mean NES x -log10FDR
  ord <- data.table(pathway = sig_any, agg = agg)
  up <- ord[agg > 0][order(-agg)]; dn <- ord[agg < 0][order(agg)]    # dn: strongest-down first
  if (is.finite(top_n_per_dir)) { up <- head(up, top_n_per_dir); dn <- head(dn, top_n_per_dir) }
  pw   <- c(up$pathway, rev(dn$pathway))                              # top = strongest up ... bottom = strongest down
  n_up <- nrow(up)

  # long table for plotting
  long <- rbindlist(lapply(comps, function(cp) {
    g <- tabs[[cp]][match(pw, pathway)]
    data.table(pathway = pw, comparison = cp, NES = g$NES, padj = g$padj)
  }))
  long[, comparison := factor(comparison, levels = comps)]
  long[, pathway    := factor(pathway, levels = rev(pw))]
  long[, star       := stars(padj)]
  long[, grey       := (is.na(padj) | padj >= sig) & (is.na(NES) | abs(NES) < grey_nes)]
  long[, NES_disp   := fifelse(grey, NA_real_, NES)]   # greyed cells render as na.value
  long[, star_col   := fifelse(!grey & !is.na(NES) & abs(NES) > 1.1, "white", "grey15")]

  fwrite(dcast(long, pathway ~ comparison, value.var = c("NES","padj")),
         file.path(out_dir, sprintf("heatmap_table_union_%s_rankBy_%s.csv", coll, rk)))

  mx <- max(abs(long$NES), na.rm = TRUE)
  p <- ggplot(long, aes(comparison, pathway, fill = NES_disp)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = star, color = star_col), size = 6, vjust = 0.78) +
    geom_vline(xintercept = burden_n + 0.5, color = "grey25", linewidth = 0.7) +
    { if (n_up > 0 && n_up < length(pw))
      geom_hline(yintercept = length(pw) - n_up + 0.5, color = "grey25", linewidth = 0.7) } +
    scale_fill_gradient2(low = "#3A5FCD", mid = "white", high = "#C44E52",
                         midpoint = 0, limits = c(-mx, mx), na.value = "grey85", name = "NES") +
    scale_color_identity() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = function(x) clean_name(x)) +
    labs(title = sprintf("Pathways FDR<%.2g in \u22651 contrast — %s | %s | rankBy %s",
                         sig, basename(pop_dir), coll, rk),
         subtitle = sprintf("rows = %d pathways (ranked by mean NES \u00d7 -log10FDR) | stars: *** <.001 ** <.01 * <.05 \u00b7 <.10 | vertical line splits burden | lenti",
                            length(pw)),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x.top = element_text(angle = 0, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))

  h <- max(3, 0.23 * length(pw) + 1.8)
  ggsave(file.path(out_dir, sprintf("heatmap_union_signature_%s_rankBy_%s.png", coll, rk)),
         p, width = 14, height = h, dpi = 150, limitsize = FALSE)
  cat(sprintf("[%s|%s] %d union-sig pathways (%d up, %d down) -> heatmap\n", coll, rk, length(pw), n_up, length(pw)-n_up))
}
cat("\nDone. Outputs in:", out_dir, "\n")
