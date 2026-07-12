## ============================================================
## H2B Label-Retaining HSC — Cross-dataset GSEA
##
## Reference: DEgene2_analysis_Res_H2BGFP_Low_VS_High.csv
##   Contrast = Low vs High
##   Positive logFC = higher in GFP-LOW  (label-LOSING / dividing HSCs)
##   Negative logFC = higher in GFP-HIGH (label-RETAINING / quiescent HSCs)
##   → ranked list is negated so positive rank = label-retaining enriched
##
## Gene sets built from ALL genes with padj < 0.05, split by direction:
##   LR_up   = padj < 0.05 AND higher in label-RETAINING HSCs (rank_score > 0)
##   LR_down = padj < 0.05 AND higher in label-LOSING HSCs    (rank_score < 0)
##
## Four LT-HSC comparisons (each tested against the H2B gene sets):
##   1. LentiNegHB_vs_LentiNegLB
##   2. LentiPosHB_vs_LentiPosLB
##   3. LentiPos_vs_LentiNeg_in_LB
##   4. LentiPos_vs_LentiNeg_in_HB
##
## OUTPUTS (written to OUTDIR):
##   GSEA_H2B_on_LentiNegHB_vs_LentiNegLB.png
##   GSEA_H2B_on_LentiPosHB_vs_LentiPosLB.png
##   GSEA_H2B_on_LentiPos_vs_LentiNeg_in_LB.png
##   GSEA_H2B_on_LentiPos_vs_LentiNeg_in_HB.png
## ============================================================

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

set.seed(42)

OUTDIR      <- "C:/Users/fc809/Downloads/H2B_fdrcutoff"
PADJ_CUTOFF <- 0.05

# Path to this script — update if you rename or move the file.
# At the end of the run, the script is copied into OUTDIR for provenance.
SCRIPT_PATH <- "C:/Users/fc809/Downloads/H2B_analysis.R"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ── FILE PATHS ────────────────────────────────────────────────────────────────

H2B_FILE <- "C:/Users/fc809/Downloads/DEgene2_analysis_Res_H2BGFP_Low_VS_High.csv"

HSC_DIR <- "C:/Users/fc809/Downloads/LT-HSCs (1)"

hsc_files <- list(
  LentiPosLB_vs_LentiNegHB   = file.path(HSC_DIR, "LentiPosLB_vs_LentiNegHB",
                                         "DEG_LentiPosLB_vs_LentiNegHB_pseudoBulk_res0.3.sce_0.txt"),
  LentiNegHB_vs_LentiNegLB   = file.path(HSC_DIR, "LentiNegHB_vs_LentiNegLB",
                                         "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  LentiPosHB_vs_LentiPosLB   = file.path(HSC_DIR, "LentiPosHB_vs_LentiPosLB",
                                         "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt"),
  LentiPos_vs_LentiNeg_in_LB = file.path(HSC_DIR, "LentiPos_vs_LentiNeg_in_LB",
                                         "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.3.sce_0.txt"),
  LentiPos_vs_LentiNeg_in_HB = file.path(HSC_DIR, "LentiPos_vs_LentiNeg_in_HB",
                                         "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.3.sce_0.txt")
)

## ── LOAD H2B DEG DATA ─────────────────────────────────────────────────────────
# Contrast is Low vs High, so positive log2FoldChange = higher in label-LOSING.
# Negate so positive rank_score = higher in label-RETAINING.

h2b <- read.csv(H2B_FILE, row.names = 1, stringsAsFactors = FALSE)
h2b$gene       <- rownames(h2b)
h2b$rank_score <- -h2b$log2FoldChange

cat("H2B DEG table loaded:", nrow(h2b), "genes\n")
cat("Rank score range:", round(range(h2b$rank_score, na.rm = TRUE), 2), "\n\n")

## ── BUILD H2B GENE SETS (all genes with padj < PADJ_CUTOFF) ───────────────────
## LR_up   = significant AND higher in label-RETAINING (rank_score > 0, i.e. log2FC < 0)
## LR_down = significant AND higher in label-LOSING    (rank_score < 0, i.e. log2FC > 0)

h2b_sig <- h2b %>%
  filter(!is.na(rank_score), !is.na(padj), padj < PADJ_CUTOFF) %>%
  distinct(gene, .keep_all = TRUE)

lr_up <- h2b_sig %>%
  filter(rank_score > 0) %>%
  arrange(desc(rank_score)) %>%
  pull(gene)

lr_down <- h2b_sig %>%
  filter(rank_score < 0) %>%
  arrange(rank_score) %>%
  pull(gene)

h2b_sets <- data.frame(
  term = c(rep("LR_up",   length(lr_up)),
           rep("LR_down", length(lr_down))),
  gene = c(lr_up, lr_down)
)

# Store set sizes for use in plot annotations later
N_LR_UP   <- length(lr_up)
N_LR_DOWN <- length(lr_down)

cat("H2B significant genes (padj <", PADJ_CUTOFF, "):", nrow(h2b_sig), "\n")
cat("  LR_up   (label-RETAINING enriched):", N_LR_UP,   "genes\n")
cat("  LR_down (label-LOSING enriched):   ", N_LR_DOWN, "genes\n\n")

## ── LOAD + RANK AN HSC DEG FILE ───────────────────────────────────────────────

detect_cols <- function(df) {
  gene_col <- colnames(df)[1]
  fc_col   <- intersect(c("logFC", "log2FoldChange", "log2FC", "LogFC"),
                        colnames(df))[1]
  if (is.na(fc_col))
    stop("Cannot find logFC column. Columns present: ",
         paste(colnames(df), collapse = ", "))
  list(gene = gene_col, fc = fc_col)
}

make_ranked_hsc <- function(filepath) {
  df   <- read.delim(filepath, stringsAsFactors = FALSE)
  cols <- detect_cols(df)
  df %>%
    filter(.data[[cols$gene]] != "Lenti") %>%
    filter(!is.na(.data[[cols$fc]])) %>%
    arrange(desc(.data[[cols$fc]])) %>%
    distinct(.data[[cols$gene]], .keep_all = TRUE) %>%
    { setNames(.[[cols$fc]], .[[cols$gene]]) }
}

## ── GSEA RUNNER ───────────────────────────────────────────────────────────────

run_gsea <- function(ranked_list) {
  # maxGSSize set high enough to retain LR_up / LR_down regardless of
  # how many genes pass the padj < PADJ_CUTOFF cutoff
  GSEA(
    geneList      = ranked_list,
    TERM2GENE     = h2b_sets,
    minGSSize     = 10,
    maxGSSize     = max(5000, N_LR_UP, N_LR_DOWN) + 1,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    verbose       = TRUE
  )
}

## ── PLOT FUNCTION ─────────────────────────────────────────────────────────────

plot_h2b_gsea <- function(gsea_obj, title, outpath) {
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  cat("\n---", title, "---\n")
  print(stats)
  
  nes_up  <- round(stats$NES[stats$ID == "LR_up"],   2)
  fdr_up  <- signif(stats$p.adjust[stats$ID == "LR_up"],   2)
  nes_dn  <- round(stats$NES[stats$ID == "LR_down"], 2)
  fdr_dn  <- signif(stats$p.adjust[stats$ID == "LR_down"], 2)
  
  p <- gseaplot2(gsea_obj,
                 geneSetID = c("LR_up", "LR_down"),
                 title     = paste0(title,
                                    "\n(H2B sets: padj<", PADJ_CUTOFF,
                                    "  |  LR_up n=", N_LR_UP,
                                    ", LR_down n=", N_LR_DOWN, ")"),
                 subplots  = 1:3,
                 base_size = 11)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(
      values = c("LR_up" = "#1A237E", "LR_down" = "#B71C1C"),
      labels = c(
        "LR_up"   = paste0("LR-up (n=", N_LR_UP,
                           ")   NES=", nes_up, "  FDR=", fdr_up),
        "LR_down" = paste0("LR-down (n=", N_LR_DOWN,
                           ") NES=", nes_dn, "  FDR=", fdr_dn)
      )
    ) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  
  p[[2]] <- p[[2]] +
    scale_color_manual(values = c("LR_up" = "#1A237E", "LR_down" = "#B71C1C"))
  
  png(outpath, width = 1400, height = 1200, res = 150)
  print(p)
  dev.off()
  message("Saved: ", basename(outpath))
}

## ── PLOT TITLES ───────────────────────────────────────────────────────────────

# The four biologically interpretable comparisons (excludes cross-group)
core_comparisons <- c(
  "LentiNegHB_vs_LentiNegLB",
  "LentiPosHB_vs_LentiPosLB",
  "LentiPos_vs_LentiNeg_in_LB",
  "LentiPos_vs_LentiNeg_in_HB"
)

plot_titles <- c(
  LentiPosLB_vs_LentiNegHB   = "H2B LR signature: LentiPos-LB vs LentiNeg-HB (LT-HSCs)",
  LentiNegHB_vs_LentiNegLB   = "H2B LR signature: LentiNeg-HB vs LentiNeg-LB (LT-HSCs)",
  LentiPosHB_vs_LentiPosLB   = "H2B LR signature: LentiPos-HB vs LentiPos-LB (LT-HSCs)",
  LentiPos_vs_LentiNeg_in_LB = "H2B LR signature: Lenti+ vs Lenti- in LB (LT-HSCs)",
  LentiPos_vs_LentiNeg_in_HB = "H2B LR signature: Lenti+ vs Lenti- in HB (LT-HSCs)"
)

# Short labels for summary plots
short_labels <- c(
  LentiNegHB_vs_LentiNegLB   = "LentiNeg\nHB vs LB",
  LentiPosHB_vs_LentiPosLB   = "LentiPos\nHB vs LB",
  LentiPos_vs_LentiNeg_in_LB = "Lenti+/-\nin LB",
  LentiPos_vs_LentiNeg_in_HB = "Lenti+/-\nin HB"
)

## ── RUN ALL FIVE — store results ──────────────────────────────────────────────

# Prevent Windows parallel worker crash in fgsea
BiocParallel::register(BiocParallel::SerialParam())

cat("=== Column detection check (first file) ===\n")
peek <- read.delim(hsc_files[[1]], nrow = 3, stringsAsFactors = FALSE)
cat("Columns:", paste(colnames(peek), collapse = ", "), "\n\n")

gsea_results  <- list()   # store GSEA objects
ranked_lists  <- list()   # store ranked gene vectors
deg_tables    <- list()   # store full DEG tables

for (nm in names(hsc_files)) {
  cat("\n=== Running:", nm, "===\n")
  deg_tables[[nm]]  <- read.delim(hsc_files[[nm]], stringsAsFactors = FALSE)
  ranked_lists[[nm]] <- make_ranked_hsc(hsc_files[[nm]])
  cat("Ranked list:", length(ranked_lists[[nm]]), "genes\n")
  gsea_results[[nm]] <- run_gsea(ranked_lists[[nm]])
  plot_h2b_gsea(
    gsea_obj = gsea_results[[nm]],
    title    = plot_titles[nm],
    outpath  = file.path(OUTDIR, paste0("GSEA_H2B_on_", nm, ".png"))
  )
}

message("\nIndividual GSEA plots done.")

## ── 1. SUMMARY DOTPLOT (four core comparisons) ────────────────────────────────

# Use clean labels without newlines for factor levels
clean_labels <- c(
  LentiNegHB_vs_LentiNegLB   = "LentiNeg HB vs LB",
  LentiPosHB_vs_LentiPosLB   = "LentiPos HB vs LB",
  LentiPos_vs_LentiNeg_in_LB = "Lenti+/- in LB",
  LentiPos_vs_LentiNeg_in_HB = "Lenti+/- in HB"
)

summary_rows <- lapply(core_comparisons, function(nm) {
  df <- as.data.frame(gsea_results[[nm]])
  df <- df[df$ID %in% c("LR_up", "LR_down"), c("ID", "NES", "p.adjust")]
  df$comparison <- clean_labels[nm]
  df
})
gsea_summary <- do.call(rbind, summary_rows)
gsea_summary$comparison <- factor(gsea_summary$comparison,
                                  levels = unname(clean_labels[core_comparisons]))
# Y-axis labels show gene-set sizes
id_labels <- c(
  LR_up   = paste0("LR_up\n(n=",   N_LR_UP,   ")"),
  LR_down = paste0("LR_down\n(n=", N_LR_DOWN, ")")
)
gsea_summary$ID <- factor(gsea_summary$ID, levels = c("LR_up", "LR_down"))

cat("\nDotplot data:\n"); print(gsea_summary)

p_dot <- ggplot(gsea_summary,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "-log10(FDR)") +
  scale_y_discrete(labels = id_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x    = element_text(size = 10, angle = 30, hjust = 1,
                                      margin = margin(t = 4)),
        axis.text.y    = element_text(size = 11),
        axis.title.y   = element_text(size = 11),
        plot.title     = element_text(size = 12, hjust = 0.5,
                                      margin = margin(b = 10)),
        plot.subtitle  = element_text(size = 10, hjust = 0.5,
                                      margin = margin(b = 8)),
        legend.margin  = margin(t = 10),
        plot.margin    = margin(t = 15, r = 10, b = 25, l = 15)) +
  labs(x = NULL, y = "Gene set",
       title    = "H2B label-retaining signature across LT-HSC comparisons",
       subtitle = paste0("Gene sets: H2B DEGs with padj < ", PADJ_CUTOFF,
                         "  |  LR_up n=", N_LR_UP,
                         ", LR_down n=", N_LR_DOWN))

ggsave(file.path(OUTDIR, "Summary_dotplot.png"),
       p_dot, width = 8, height = 4.5, dpi = 150)
message("Saved: Summary_dotplot.png")

## ── 2. LEADING EDGE HEATMAP ───────────────────────────────────────────────────

library(pheatmap)

# Robust leading edge extraction — handles both character and list-column storage
extract_leading_edge <- function(gsea_obj, set_id) {
  df <- as.data.frame(gsea_obj)
  row <- df[df$ID == set_id, , drop = FALSE]
  if (nrow(row) == 0) return(character(0))
  le <- row$core_enrichment[1]
  if (is.na(le) || le == "") return(character(0))
  trimws(unlist(strsplit(as.character(le), "/")))
}

# Diagnostic: show leading edge sizes
cat("\n=== Leading edge sizes ===\n")
for (nm in core_comparisons) {
  for (sid in c("LR_up", "LR_down")) {
    le <- extract_leading_edge(gsea_results[[nm]], sid)
    cat(nm, "|", sid, ":", length(le), "genes\n")
  }
}

cols_info <- detect_cols(deg_tables[[core_comparisons[1]]])
gene_col  <- cols_info$gene
fc_col    <- cols_info$fc

for (set_id in c("LR_up", "LR_down")) {
  
  le_genes <- unique(unlist(lapply(core_comparisons, function(nm)
    extract_leading_edge(gsea_results[[nm]], set_id))))
  
  cat("\nLeading edge union (", set_id, "):", length(le_genes), "genes\n")
  
  if (length(le_genes) == 0) {
    message("No leading edge genes found for ", set_id, " — skipping heatmap")
    next
  }
  
  # Build logFC matrix (genes x comparisons)
  fc_mat <- sapply(core_comparisons, function(nm) {
    df <- deg_tables[[nm]]
    fc <- setNames(df[[fc_col]], df[[gene_col]])
    fc[le_genes]
  })
  rownames(fc_mat) <- le_genes
  colnames(fc_mat) <- unname(clean_labels[core_comparisons])
  
  # Drop genes missing in more than half the comparisons
  fc_mat <- fc_mat[rowSums(is.na(fc_mat)) <= 2, , drop = FALSE]
  fc_mat[is.na(fc_mat)] <- 0
  # Cap extreme values for color scale
  fc_mat <- pmax(pmin(fc_mat, 2), -2)
  
  cat("Heatmap matrix:", nrow(fc_mat), "x", ncol(fc_mat), "\n")
  
  png(file.path(OUTDIR, paste0("LeadingEdge_heatmap_", set_id, ".png")),
      width = 1000, height = max(700, nrow(fc_mat) * 14 + 250), res = 130)
  set_total <- if (set_id == "LR_up") N_LR_UP else N_LR_DOWN
  pheatmap(fc_mat,
           color             = colorRampPalette(c("#1A237E", "white", "#B71C1C"))(100),
           breaks            = seq(-2, 2, length.out = 101),
           cluster_cols      = FALSE,
           clustering_method = "ward.D2",
           fontsize_row      = 7,
           fontsize_col      = 12,
           border_color      = NA,
           main              = paste0("Leading edge genes — ", set_id,
                                      "  (", nrow(fc_mat), " / ", set_total,
                                      " set genes)"))
  dev.off()
  message("Saved: LeadingEdge_heatmap_", set_id, ".png")
}

## ── 3. LEADING EDGE VENN DIAGRAMS ─────────────────────────────────────────────

library(ggVennDiagram)

for (set_id in c("LR_up", "LR_down")) {
  
  le_list <- setNames(
    lapply(core_comparisons, function(nm)
      extract_leading_edge(gsea_results[[nm]], set_id)),
    unname(clean_labels[core_comparisons])
  )
  
  # Skip if all empty
  if (all(lengths(le_list) == 0)) {
    message("No leading edge genes for ", set_id, " — skipping Venn")
    next
  }
  
  set_total <- if (set_id == "LR_up") N_LR_UP else N_LR_DOWN
  p_venn <- ggVennDiagram(le_list, label_alpha = 0, set_size = 3.5) +
    scale_fill_gradient(low = "white",
                        high = ifelse(set_id == "LR_up", "#1A237E", "#B71C1C")) +
    labs(title = paste0("Leading edge overlap — ", set_id,
                        "  (of ", set_total, " set genes)")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 13))
  
  ggsave(file.path(OUTDIR, paste0("LeadingEdge_Venn_", set_id, ".png")),
         p_venn, width = 7, height = 6, dpi = 150)
  message("Saved: LeadingEdge_Venn_", set_id, ".png")
}

## ── 4. RANKED LIST CORRELATION SCATTERS ───────────────────────────────────────

scatter_pairs <- list(
  list(x = "LentiNegHB_vs_LentiNegLB",
       y = "LentiPosHB_vs_LentiPosLB",
       xlab = "logFC (LentiNeg HB vs LB)",
       ylab = "logFC (LentiPos HB vs LB)",
       title = "HB vs LB: LentiNeg vs LentiPos"),
  list(x = "LentiPos_vs_LentiNeg_in_LB",
       y = "LentiPos_vs_LentiNeg_in_HB",
       xlab = "logFC (Lenti+/- in LB)",
       ylab = "logFC (Lenti+/- in HB)",
       title = "Lenti+/- effect: LB vs HB context")
)

all_le <- unique(unlist(lapply(core_comparisons, function(nm)
  c(extract_leading_edge(gsea_results[[nm]], "LR_up"),
    extract_leading_edge(gsea_results[[nm]], "LR_down")))))
cat("\nTotal leading edge genes for scatter highlight:", length(all_le), "\n")

plot_list <- lapply(scatter_pairs, function(p) {
  fc_x   <- ranked_lists[[p$x]]
  fc_y   <- ranked_lists[[p$y]]
  shared <- intersect(names(fc_x), names(fc_y))
  
  df <- data.frame(
    gene  = shared,
    x     = as.numeric(fc_x[shared]),
    y     = as.numeric(fc_y[shared]),
    is_le = shared %in% all_le
  )
  df <- df[complete.cases(df), ]
  
  r <- sprintf("%.3f", cor(df$x, df$y))
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data = subset(df, !is_le),
               color = "grey80", size = 0.4, alpha = 0.5) +
    geom_point(data = subset(df,  is_le),
               aes(color = x + y), size = 1.5, alpha = 0.9) +
    scale_color_gradient2(low = "#1A237E", mid = "grey50", high = "#B71C1C",
                          midpoint = 0,
                          name = "Leading edge
genes
(sum of logFC)",
                          breaks = c(-2, 0, 2),
                          labels = c("Down in both", "Mixed", "Up in both")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_smooth(method = "lm", se = FALSE, color = "black",
                linewidth = 0.6, linetype = "dotted") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", r), size = 4) +
    labs(x = p$xlab, y = p$ylab, title = p$title) +
    guides(color = guide_colorbar(barwidth = 0.8, barheight = 4,
                                  title.position = "top")) +
    theme_classic(base_size = 11) +
    theme(legend.title    = element_text(size = 8),
          legend.text     = element_text(size = 7),
          legend.position = c(0.13, 0.82),
          legend.background = element_rect(fill = alpha("white", 0.7),
                                           color = NA))
})

p_scatter <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
ggsave(file.path(OUTDIR, "RankedList_scatter.png"),
       p_scatter, width = 11, height = 5, dpi = 150)
message("Saved: RankedList_scatter.png")

message("\n=== All outputs saved to: ", OUTDIR, " ===")
message("  Individual GSEA curves:  5 PNGs")
message("  Summary dotplot:         Summary_dotplot.png")
message("  Leading edge heatmaps:   LeadingEdge_heatmap_LR_up/down.png")
message("  Leading edge Venns:      LeadingEdge_Venn_LR_up/down.png")
message("  Ranked list scatters:    RankedList_scatter.png")

## ============================================================
## H2B LABEL-RETAINING vs LABEL-LOSING HSC — PATHWAY ANALYSIS
##
## Uses the full H2B DEG ranked list (negated: positive = LR-enriched)
## Three analyses:
##   A. GSEA — GO Biological Process (gseGO)
##   B. GSEA — KEGG (gseKEGG)
##   C. GSEA — MSigDB Hallmarks
##
## OUTPUTS (written to OUTDIR/H2B_pathways/):
##   H2B_GSEA_GOBP_bar.png
##   H2B_GSEA_GOBP_dot.png
##   H2B_GSEA_KEGG_bar.png
##   H2B_GSEA_KEGG_dot.png
##   H2B_GSEA_Hallmarks_bar.png
##   H2B_GSEA_Hallmarks_dot.png
##   H2B_GSEA_GOBP.csv
##   H2B_GSEA_KEGG.csv
##   H2B_GSEA_Hallmarks.csv
## ============================================================

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)

PW_DIR <- file.path(OUTDIR, "H2B_pathways")
dir.create(PW_DIR, showWarnings = FALSE)

## ── Ranked list: positive = higher in label-RETAINING ────────────────────────
# h2b$rank_score already computed above (-log2FoldChange)

ranked_h2b_full <- h2b %>%
  filter(!is.na(rank_score)) %>%
  arrange(desc(rank_score)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$rank_score, .$gene) }

cat("Full H2B ranked list:", length(ranked_h2b_full), "genes\n")

## ── Convert gene symbols to Entrez IDs (required for KEGG) ───────────────────
sym2entrez <- bitr(names(ranked_h2b_full),
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Mm.eg.db)
cat("Mapped", nrow(sym2entrez), "of", length(ranked_h2b_full), "genes to Entrez\n")

ranked_h2b_entrez <- ranked_h2b_full[sym2entrez$SYMBOL]
names(ranked_h2b_entrez) <- sym2entrez$ENTREZID
# Keep highest |score| if duplicate Entrez IDs
ranked_h2b_entrez <- ranked_h2b_entrez[!duplicated(names(ranked_h2b_entrez))]
ranked_h2b_entrez <- sort(ranked_h2b_entrez, decreasing = TRUE)

## ── Plot helper ───────────────────────────────────────────────────────────────

save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
  
  res_df <- as.data.frame(gsea_obj)
  write.csv(res_df, file.path(PW_DIR, paste0(prefix, ".csv")), row.names = FALSE)
  cat("\n---", prefix, "--- significant (FDR<0.05):", sum(res_df$p.adjust < 0.05), "\n")
  
  sig <- res_df[res_df$p.adjust < 0.05, ]
  if (nrow(sig) == 0) {
    message("No significant pathways for ", prefix, " at FDR<0.05 — using top ", top_n)
    sig <- res_df[order(res_df$pvalue), ][seq_len(min(top_n, nrow(res_df))), ]
  }
  
  # Barplot: top enriched and depleted
  top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
  top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
  bar_df  <- rbind(top_pos, top_neg)
  
  if (nrow(bar_df) > 0) {
    # Trim long descriptions
    bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                 paste0(substr(bar_df$Description, 1, 52), "..."),
                                 bar_df$Description)
    bar_df$Description <- factor(bar_df$Description,
                                 levels = bar_df$Description[order(bar_df$NES)])
    
    p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                          limits = c(0, 0.05)) +
      geom_vline(xintercept = 0, color = "grey30") +
      labs(x = "Normalised Enrichment Score",
           y = NULL,
           title = paste0(prefix, "  |  positive NES = label-retaining enriched")) +
      theme_classic(base_size = 10) +
      theme(plot.title  = element_text(size = 9, hjust = 0.5),
            axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PW_DIR, paste0(prefix, "_bar.png")),
           p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5), dpi = 150)
    message("Saved: ", prefix, "_bar.png")
  }
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      theme(strip.text  = element_text(size = 9),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(PW_DIR, paste0(prefix, "_dot.png")),
           p_dot, width = 10,
           height = max(5, min(top_n, nrow(res_df)) * 0.28 + 2), dpi = 150)
    message("Saved: ", prefix, "_dot.png")
  }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
}

## ── A. GO Biological Process ─────────────────────────────────────────────────

cat("\n=== GO BP ===\n")
gsea_go <- gseGO(
  geneList     = ranked_h2b_full,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps          = 0,
  verbose      = TRUE
)
gsea_go <- setReadable(gsea_go, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
save_pathway_plots(gsea_go, "H2B_GSEA_GOBP")

## ── B. KEGG ───────────────────────────────────────────────────────────────────

cat("\n=== KEGG ===\n")
gsea_kegg <- gseKEGG(
  geneList      = ranked_h2b_entrez,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save_pathway_plots(gsea_kegg, "H2B_GSEA_KEGG")

## ── C. MSigDB Hallmarks ───────────────────────────────────────────────────────

cat("\n=== Hallmarks ===\n")
hallmarks <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks <- GSEA(
  geneList      = ranked_h2b_full,
  TERM2GENE     = hallmarks,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
save_pathway_plots(gsea_hallmarks, "H2B_GSEA_Hallmarks")

message("\n=== Pathway analysis complete. Results in: ", PW_DIR, " ===")

## ── COPY THIS SCRIPT INTO OUTDIR FOR PROVENANCE ───────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  script_dest <- file.path(OUTDIR, basename(SCRIPT_PATH))
  file.copy(SCRIPT_PATH, script_dest, overwrite = TRUE)
  message("Script copied to: ", script_dest)
} else {
  warning("SCRIPT_PATH not found — script not copied into OUTDIR: ", SCRIPT_PATH)
}

