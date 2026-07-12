# ══════════════════════════════════════════════════════════════════════════════
#  Summary Bar Plots
#    1. DEGs (FDR < 0.2): excl. HB3/HB4 vs all samples
#    2. Significant pathways (FDR < 0.05) per comparison x database
#       for both all samples and excl. HB3/HB4
#  For both LT-HSCs and Classical Monocytes
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

# ── Comparisons ───────────────────────────────────────────────────────────────
comparisons <- list(
  list(folder = "LentiPos_vs_LentiNeg_in_HB", label = "LentiPos vs LentiNeg\n(HB)"),
  list(folder = "LentiPos_vs_LentiNeg_in_LB", label = "LentiPos vs LentiNeg\n(LB)"),
  list(folder = "LentiPosHB_vs_LentiPosLB",   label = "LentiPos\nHB vs LB"),
  list(folder = "LentiNegHB_vs_LentiNegLB",   label = "LentiNeg\nHB vs LB")
)
comp_labels <- map_chr(comparisons, "label")

# ── Databases ─────────────────────────────────────────────────────────────────
db_tags <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Hallmark")
db_colors <- c(
  GO_BP    = "#3498DB",
  GO_MF    = "#2ECC71",
  GO_CC    = "#1ABC9C",
  KEGG     = "#E74C3C",
  Hallmark = "#9B59B6"
)

# ── Datasets ──────────────────────────────────────────────────────────────────
datasets <- list(
  list(
    name     = "LT-HSCs",
    all_dir  = "LT-HSCs (1)/allSamples_analysis",
    excl_dir = "LT-HSCs (1)/excluding_H3_H4",
    out      = "LT-HSCs (1)/summary_plots"
  ),
  list(
    name     = "Classical Monocytes",
    all_dir  = "Classical_Monocytes (1)/allSamples_analysis",
    excl_dir = "Classical_Monocytes (1)/excluding_H3_H4",
    out      = "Classical_Monocytes (1)/summary_plots"
  )
)

DEG_FDR      <- 0.1
PATHWAY_FDR  <- 0.05

# ── Helpers ───────────────────────────────────────────────────────────────────
find_limma <- function(base, folder) {
  p     <- file.path(base, folder)
  files <- list.files(p, pattern = "_limma_.*results\\.csv$", full.names = TRUE)
  if (length(files) > 0) return(files[1])
  NULL
}

find_gsea <- function(base, folder, db) {
  p     <- file.path(base, folder)
  files <- list.files(p, pattern = paste0("_GSEA_", db, "\\.csv$"), full.names = TRUE)
  if (length(files) > 0) return(files[1])
  NULL
}

load_limma <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  if (!"adj.P.Val" %in% names(df) && "FDR"  %in% names(df)) df <- dplyr::rename(df, adj.P.Val = "FDR")
  if (!"adj.P.Val" %in% names(df) && "padj" %in% names(df)) df <- dplyr::rename(df, adj.P.Val = "padj")
  df
}

# ══════════════════════════════════════════════════════════════════════════════
for (ds in datasets) {
  dir.create(ds$out, showWarnings = FALSE, recursive = TRUE)
  cat("\n══", ds$name, "══\n")
  
  analyses <- list(
    list(label = "All samples",   dir = ds$all_dir),
    list(label = "Excl. HB3/HB4", dir = ds$excl_dir)
  )
  
  # ── 1. DEG counts ──────────────────────────────────────────────────────────
  deg_rows <- list()
  for (an in analyses) {
    for (comp in comparisons) {
      f <- find_limma(an$dir, comp$folder)
      if (is.null(f)) { cat("  [DEG] not found:", an$label, "|", comp$folder, "\n"); next }
      res    <- load_limma(f)
      n_up   <- sum(res$adj.P.Val < DEG_FDR & res$logFC >  0, na.rm = TRUE)
      n_down <- sum(res$adj.P.Val < DEG_FDR & res$logFC <  0, na.rm = TRUE)
      deg_rows[[paste(an$label, comp$folder)]] <- tibble(
        analysis   = an$label,
        comparison = comp$label,
        direction  = c("Up", "Down"),
        n_genes    = c(n_up, n_down)
      )
      cat("  [DEG]", an$label, "|", comp$folder, "— up:", n_up, "down:", n_down, "\n")
    }
  }
  
  if (length(deg_rows) > 0) {
    deg_df <- bind_rows(deg_rows) %>%
      mutate(
        analysis   = factor(analysis,   levels = c("All samples", "Excl. HB3/HB4")),
        comparison = factor(comparison, levels = comp_labels),
        direction  = factor(direction,  levels = c("Up", "Down"))
      )
    
    p_deg <- ggplot(deg_df, aes(x = comparison, y = n_genes, fill = direction)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
      geom_text(aes(label = n_genes),
                position = position_dodge(width = 0.7),
                vjust = -0.4, size = 3.5, fontface = "bold", color = "black") +
      facet_wrap(~ analysis, nrow = 2) +
      scale_fill_manual(values = c(Up = "#E74C3C", Down = "#3498DB"), name = "Direction") +
      labs(
        title    = paste0(ds$name, " — Differentially Expressed Genes (FDR < ", DEG_FDR, ")"),
        subtitle = "Top: all samples  |  Bottom: excluding HB3/HB4",
        x = NULL, y = "Number of DEGs"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title         = element_text(face = "bold", size = 14, color = "black"),
        plot.subtitle      = element_text(size = 11, color = "grey30"),
        axis.text.x        = element_text(size = 11, color = "black", lineheight = 1.2),
        axis.text.y        = element_text(size = 11, color = "black"),
        axis.title.y       = element_text(size = 12, color = "black"),
        strip.text         = element_text(size = 11, face = "bold"),
        strip.background   = element_rect(fill = "grey92", color = NA),
        legend.text        = element_text(size = 11),
        legend.title       = element_text(size = 11, face = "bold"),
        panel.grid.major.x = element_blank()
      )
    
    ggsave(file.path(ds$out, "DEG_counts_FDR0.1.pdf"), p_deg, width = 13, height = 9)
    ggsave(file.path(ds$out, "DEG_counts_FDR0.1.png"), p_deg, width = 13, height = 9, dpi = 160)
    cat("  Saved: DEG_counts_FDR0.1\n")
  }
  
  # ── 2. Pathway counts — both analyses ─────────────────────────────────────
  for (an in analyses) {
    pathway_rows <- list()
    for (comp in comparisons) {
      for (db in db_tags) {
        f <- find_gsea(an$dir, comp$folder, db)
        if (is.null(f)) {
          pathway_rows[[paste(comp$folder, db)]] <- tibble(
            comparison = comp$label, database = db,
            n_sig = 0L, n_up = 0L, n_down = 0L
          )
          next
        }
        res <- read_csv(f, show_col_types = FALSE)
        if (!"p.adjust" %in% names(res) && "qvalue" %in% names(res))
          res <- dplyr::rename(res, p.adjust = "qvalue")
        pathway_rows[[paste(comp$folder, db)]] <- tibble(
          comparison = comp$label, database = db,
          n_sig  = sum(res$p.adjust < PATHWAY_FDR, na.rm = TRUE),
          n_up   = sum(res$p.adjust < PATHWAY_FDR & res$NES > 0, na.rm = TRUE),
          n_down = sum(res$p.adjust < PATHWAY_FDR & res$NES < 0, na.rm = TRUE)
        )
      }
    }
    
    pathway_df <- bind_rows(pathway_rows) %>%
      mutate(
        comparison = factor(comparison, levels = comp_labels),
        database   = factor(database,   levels = db_tags)
      )
    
    an_tag <- gsub("[^A-Za-z0-9]", "_", an$label)
    
    # Stacked bar
    p_stack <- ggplot(pathway_df, aes(x = comparison, y = n_sig, fill = database)) +
      geom_col(width = 0.65, alpha = 0.85) +
      geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")),
                position = position_stack(vjust = 0.5),
                size = 3.2, fontface = "bold", color = "white") +
      scale_fill_manual(values = db_colors, name = "Database") +
      labs(
        title    = paste0(ds$name, " — Significant Pathways (FDR < ", PATHWAY_FDR, ")"),
        subtitle = paste0(an$label, "  |  Stacked by database"),
        x = NULL, y = "Number of significant pathways"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title         = element_text(face = "bold", size = 14, color = "black"),
        plot.subtitle      = element_text(size = 11, color = "grey30"),
        axis.text.x        = element_text(size = 11, color = "black", lineheight = 1.2),
        axis.text.y        = element_text(size = 11, color = "black"),
        axis.title.y       = element_text(size = 12, color = "black"),
        legend.text        = element_text(size = 11),
        legend.title       = element_text(size = 11, face = "bold"),
        panel.grid.major.x = element_blank()
      )
    
    # Faceted by database, split up/down
    pathway_long <- pathway_df %>%
      dplyr::select(comparison, database, n_up, n_down) %>%
      pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n") %>%
      mutate(
        direction = recode(direction, n_up = "Up (NES > 0)", n_down = "Down (NES < 0)"),
        direction = factor(direction, levels = c("Up (NES > 0)", "Down (NES < 0)"))
      )
    
    p_facet <- ggplot(pathway_long, aes(x = comparison, y = n, fill = direction)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
      geom_text(aes(label = ifelse(n > 0, n, "")),
                position = position_dodge(width = 0.7),
                vjust = -0.4, size = 3.2, fontface = "bold", color = "black") +
      facet_wrap(~ database, scales = "free_y", nrow = 1) +
      scale_fill_manual(
        values = c("Up (NES > 0)" = "#E74C3C", "Down (NES < 0)" = "#3498DB"),
        name = "Direction"
      ) +
      labs(
        title    = paste0(ds$name, " — Significant Pathways by Database (FDR < ", PATHWAY_FDR, ")"),
        subtitle = paste0(an$label, "  |  Split by enrichment direction"),
        x = NULL, y = "Number of significant pathways"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title         = element_text(face = "bold", size = 14, color = "black"),
        plot.subtitle      = element_text(size = 11, color = "grey30"),
        axis.text.x        = element_text(size = 9, color = "black", angle = 35,
                                          hjust = 1, lineheight = 1.1),
        axis.text.y        = element_text(size = 10, color = "black"),
        axis.title.y       = element_text(size = 11, color = "black"),
        strip.text         = element_text(size = 11, face = "bold"),
        strip.background   = element_rect(fill = "grey92", color = NA),
        legend.text        = element_text(size = 10),
        legend.title       = element_text(size = 10, face = "bold"),
        panel.grid.major.x = element_blank()
      )
    
    ggsave(file.path(ds$out, paste0("pathway_counts_stacked_",     an_tag, ".pdf")), p_stack, width = 13, height = 6)
    ggsave(file.path(ds$out, paste0("pathway_counts_stacked_",     an_tag, ".png")), p_stack, width = 13, height = 6, dpi = 160)
    ggsave(file.path(ds$out, paste0("pathway_counts_by_database_", an_tag, ".pdf")), p_facet, width = 18, height = 6)
    ggsave(file.path(ds$out, paste0("pathway_counts_by_database_", an_tag, ".png")), p_facet, width = 18, height = 6, dpi = 160)
    cat("  Saved: pathway_counts —", an$label, "\n")
  }
}

cat("\n✓ Done.\n")

