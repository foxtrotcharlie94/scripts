###artifact_analysis####
#1_PCA#
# ============================================================
# PCA — Pseudobulk RNA-seq (CPM input)
# HB vs LB_LentiPos / LentiNeg_LT-HSCs & Classical Monocytes
# ============================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
setwd("C:/Users/fc809/Downloads")

files <- list(
  "LT-HSC_LentiNeg"    = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSC_LentiPos"    = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "Monocytes_LentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt",
  "Monocytes_LentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt"
)

OUTLIERS <- c("HB3", "HB4")

# ── Helper: run PCA for one file ─────────────────────────────
# Input is CPM: apply log2(CPM+1) before PCA
run_pca <- function(path, dataset_name) {
  df <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  count_cols <- grep("HB|LB", colnames(df), value = TRUE)
  cpm        <- df[, count_cols, drop = FALSE]   # already CPM
  
  # log2(CPM + 1), then z-score genes; samples as observations
  mat        <- t(log2(cpm + 1))
  mat_scaled <- scale(mat)
  
  pca     <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  var_exp <- summary(pca)$importance[2, 1:2]
  
  coords         <- as.data.frame(pca$x[, 1:2])
  coords$sample  <- rownames(coords)
  coords$label   <- sub(".*(HB|LB)(\\d)", "\\1\\2", coords$sample)
  coords$exposure<- ifelse(grepl("HB", coords$sample), "HB", "LB")
  coords$outlier <- coords$label %in% OUTLIERS & coords$exposure == "HB"
  coords$dataset <- dataset_name
  coords$var_pc1 <- var_exp[1]
  coords$var_pc2 <- var_exp[2]
  coords
}

pca_list <- mapply(run_pca, files, names(files), SIMPLIFY = FALSE)

make_pca_plot <- function(df) {
  xlab <- sprintf("PC1 (%.1f%%)", df$var_pc1[1] * 100)
  ylab <- sprintf("PC2 (%.1f%%)", df$var_pc2[1] * 100)
  
  df$color_group <- ifelse(df$exposure == "LB", "LB",
                           ifelse(df$outlier,           "HB (outlier)", "HB"))
  
  ggplot(df, aes(x = PC1, y = PC2, color = color_group, label = label)) +
    geom_point(data = subset(df, outlier), aes(x = PC1, y = PC2),
               shape = 21, size = 10, stroke = 1.2,
               color = "#ff5252", fill = NA, linetype = "dashed", inherit.aes = FALSE) +
    geom_point(aes(shape = exposure), size = 4.5, alpha = 0.9) +
    geom_text_repel(size = 3, fontface = "bold", min.segment.length = 0,
                    box.padding = 0.4, show.legend = FALSE) +
    scale_color_manual(name = NULL,
                       values = c("HB" = "#4fc3f7", "HB (outlier)" = "#ff5252", "LB" = "#69db7c")) +
    scale_shape_manual(name = NULL, values = c("HB" = 19, "LB" = 17)) +
    labs(title = unique(df$dataset), x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.key.size = unit(0.5, "cm"))
}

plots    <- lapply(pca_list, make_pca_plot)
combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "Pseudobulk PCA — HB vs LB across compartments",
                  theme = theme(plot.title = element_text(face = "bold", size = 13),
                                plot.subtitle = element_text(size = 10, color = "firebrick")))

ggsave("pca_pseudobulk.pdf", combined, width = 12, height = 10)
ggsave("pca_pseudobulk.png", combined, width = 12, height = 10, dpi = 300)
message("Saved: pca_pseudobulk.pdf / pca_pseudobulk.png")

# ── Optional: inspect top genes driving PC1 ───────────────────
inspect_loadings <- function(path, n_genes = 20) {
  df         <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  count_cols <- grep("HB|LB", colnames(df), value = TRUE)
  mat        <- t(log2(df[, count_cols] + 1))
  mat_scaled <- scale(mat)
  pca        <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  loadings   <- sort(pca$rotation[, 1], decreasing = TRUE)
  cat("\nTop genes with HIGHEST PC1 loading:\n"); print(head(loadings, n_genes))
  cat("\nTop genes with LOWEST PC1 loading:\n");  print(tail(loadings, n_genes))
}
# Uncomment to run:
# inspect_loadings(files[["LT-HSC_LentiNeg"]])


########GSEA_LIMMA#########
#PATHWAYS

# ============================================================
# GSEA Volcano Plots — HB3+HB4 vs HB1+HB2+HB5  (CPM input)
# All four datasets: LT-HSC & Monocytes × LentiNeg & LentiPos
# DEG via limma on log2(CPM+1) — no voom, no TMM
# ============================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
library(limma)
library(fgsea)
library(msigdbr)

files <- list(
  "LT-HSC\nLentiNeg"    = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSC\nLentiPos"    = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "Monocytes\nLentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt",
  "Monocytes\nLentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt"
)

message("Loading MSigDB Hallmark gene sets (mouse)...")
msig     <- msigdbr(species = "Mus musculus", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

clean_name <- function(x) {
  x <- sub("^HALLMARK_", "", x); x <- gsub("_", " ", x); tools::toTitleCase(tolower(x))
}

# ── DEG: log2(CPM+1), limma without voom ─────────────────────
run_deg <- function(path, dataset_name) {
  message("\nProcessing: ", dataset_name)
  
  df          <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  hb_cols     <- grep("HB", colnames(df), value = TRUE)
  cpm_mat     <- df[, hb_cols, drop = FALSE]
  
  sample_nums <- sub(".*HB(\\d)$", "\\1", hb_cols)
  group       <- factor(ifelse(sample_nums %in% c("3", "4"), "outlier", "normal"),
                        levels = c("normal", "outlier"))
  
  message("  Samples: ", paste(hb_cols, collapse = ", "))
  message("  Groups:  ", paste(group,   collapse = ", "))
  
  # Filter: keep genes with mean CPM >= 1 in at least one group
  grp_means <- sapply(levels(group), function(g)
    rowMeans(cpm_mat[, group == g, drop = FALSE]))
  keep    <- rowSums(grp_means >= 1) >= 1
  cpm_mat <- cpm_mat[keep, ]
  message("  Genes after filtering: ", nrow(cpm_mat))
  
  log_cpm <- log2(cpm_mat + 1)
  design  <- model.matrix(~ group)
  fit     <- eBayes(lmFit(log_cpm, design))
  
  res      <- topTable(fit, coef = "groupoutlier", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res
}

# ── GSEA ─────────────────────────────────────────────────────
run_gsea <- function(deg_res) {
  ranks <- setNames(sign(deg_res$logFC) * -log10(deg_res$P.Value), deg_res$gene)
  ranks <- sort(ranks[!is.na(ranks) & is.finite(ranks)], decreasing = TRUE)
  set.seed(42)
  fgsea(pathways = pathways, stats = ranks,
        minSize = 10, maxSize = 500, nPermSimple = 10000)
}

# ── Volcano plot ─────────────────────────────────────────────
make_gsea_volcano <- function(gsea_res, title, n_label = 12) {
  df               <- as.data.frame(gsea_res)
  df$pathway_clean <- clean_name(df$pathway)
  df$log10_padj    <- -log10(df$padj + 1e-10)
  df$sig           <- df$padj < 0.05
  df$pt_color      <- ifelse(!df$sig, "ns", ifelse(df$NES > 0, "up", "down"))
  
  up   <- df[df$sig & df$NES > 0, ]; up   <- up[order(up$padj),   ][seq_len(min(ceiling(n_label/2), nrow(up))),   ]
  down <- df[df$sig & df$NES < 0, ]; down <- down[order(down$padj),][seq_len(min(floor(n_label/2),  nrow(down))), ]
  label_df <- if (nrow(up) + nrow(down) < 4) df[order(df$padj), ][seq_len(min(n_label, nrow(df))), ] else rbind(up, down)
  
  col_vals <- c(up = "#c0392b", down = "#2980b9", ns = "#bdc3c7")
  
  ggplot(df, aes(x = NES, y = log10_padj)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "#ccc", linewidth = 0.5) +
    geom_point(aes(color = pt_color, size = log10_padj), alpha = 0.8) +
    geom_text_repel(data = label_df, aes(x = NES, y = log10_padj, label = pathway_clean, color = pt_color),
                    inherit.aes = FALSE, size = 2.6, fontface = "bold",
                    max.overlaps = 25, box.padding = 0.45, min.segment.length = 0.2, show.legend = FALSE) +
    scale_color_manual(values = col_vals,
                       labels = c(up = "Enriched in HB3/4", down = "Depleted in HB3/4", ns = "n.s."), name = NULL) +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    labs(title = title, x = "NES  (positive = enriched in HB3/4)", y = expression(-log[10](padj))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 9))
}

deg_results  <- mapply(run_deg,  files, names(files), SIMPLIFY = FALSE)
gsea_results <- lapply(deg_results, run_gsea)
plots        <- mapply(make_gsea_volcano, gsea_results, names(files), SIMPLIFY = FALSE)

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "GSEA Volcano - HB3+HB4 vs HB1+HB2+HB5",
                  subtitle = "Hallmark gene sets (MSigDB) | Ranked by sign(logFC) x -log10(P) | limma on log2(CPM+1)",
                  theme = theme(plot.title    = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 10, color = "grey40")))

ggsave("gsea_volcano_HB34_vs_HB125.pdf", combined, width = 14, height = 11)
ggsave("gsea_volcano_HB34_vs_HB125.png", combined, width = 14, height = 11, dpi = 300)
message("Saved!")

for (nm in names(gsea_results)) {
  cat("\n===", nm, "===\n")
  top <- as.data.frame(gsea_results[[nm]])
  top <- top[order(top$padj), ][1:10, c("pathway", "NES", "padj")]
  top$pathway <- clean_name(top$pathway)
  print(top, row.names = FALSE)
}


#INDIVIDUAL_GENES
# ============================================================
# Gene Volcano Plots — HB3+HB4 vs HB1+HB2+HB5  (CPM input)
# limma on log2(CPM+1), no voom, no TMM
# ============================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
library(limma)

files <- list(
  "LT-HSC\nLentiNeg"    = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSC\nLentiPos"    = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "Monocytes\nLentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt",
  "Monocytes\nLentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt"
)

FC_THRESH  <- 0.5
FDR_THRESH <- 0.05
N_LABEL    <- 20

run_deg <- function(path, dataset_name) {
  message("\nProcessing: ", dataset_name)
  
  df          <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  hb_cols     <- grep("HB", colnames(df), value = TRUE)
  cpm_mat     <- df[, hb_cols, drop = FALSE]
  
  sample_nums <- sub(".*HB(\\d)$", "\\1", hb_cols)
  group       <- factor(ifelse(sample_nums %in% c("3", "4"), "outlier", "normal"),
                        levels = c("normal", "outlier"))
  
  message("  Samples: ", paste(hb_cols, collapse = ", "))
  message("  Groups:  ", paste(group,   collapse = ", "))
  
  grp_means <- sapply(levels(group), function(g)
    rowMeans(cpm_mat[, group == g, drop = FALSE]))
  keep    <- rowSums(grp_means >= 1) >= 1
  cpm_mat <- cpm_mat[keep, ]
  message("  Genes after filtering: ", nrow(cpm_mat))
  
  log_cpm <- log2(cpm_mat + 1)
  design  <- model.matrix(~ group)
  fit     <- eBayes(lmFit(log_cpm, design))
  
  res      <- topTable(fit, coef = "groupoutlier", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res
}

make_gene_volcano <- function(deg_res, title,
                              fc_thresh  = FC_THRESH,
                              fdr_thresh = FDR_THRESH,
                              n_label    = N_LABEL) {
  df           <- deg_res
  df$log10_fdr <- -log10(df$adj.P.Val + 1e-300)
  df$sig       <- ifelse(df$adj.P.Val < fdr_thresh & df$logFC >  fc_thresh, "up",
                         ifelse(df$adj.P.Val < fdr_thresh & df$logFC < -fc_thresh, "down", "ns"))
  
  n_up   <- sum(df$sig == "up")
  n_down <- sum(df$sig == "down")
  n_ns   <- sum(df$sig == "ns")
  
  up_lab   <- df[df$sig == "up",   ]; up_lab   <- up_lab[order(up_lab$adj.P.Val),   ][seq_len(min(ceiling(n_label/2), nrow(up_lab))),   ]
  down_lab <- df[df$sig == "down", ]; down_lab <- down_lab[order(down_lab$adj.P.Val),][seq_len(min(floor(n_label/2),  nrow(down_lab))), ]
  label_df <- if (nrow(up_lab) + nrow(down_lab) == 0) df[order(df$adj.P.Val), ][seq_len(min(n_label, nrow(df))), ] else rbind(up_lab, down_lab)
  
  col_vals <- c(up = "#c0392b", down = "#2980b9", ns = "#ccc")
  subtitle <- sprintf("Up in HB3/4: %d    Down in HB3/4: %d    n.s.: %d  |  FDR < %.2f, |logFC| > %.1f",
                      n_up, n_down, n_ns, fdr_thresh, fc_thresh)
  
  ggplot(df, aes(x = logFC, y = log10_fdr)) +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "#999", linewidth = 0.5) +
    geom_vline(xintercept =  fc_thresh, linetype = "dashed", color = "#999", linewidth = 0.5) +
    geom_vline(xintercept = -fc_thresh, linetype = "dashed", color = "#999", linewidth = 0.5) +
    geom_point(data = df[df$sig == "ns",  ], aes(color = sig), size = 1,   alpha = 0.4) +
    geom_point(data = df[df$sig != "ns",  ], aes(color = sig), size = 1.8, alpha = 0.85) +
    geom_text_repel(data = label_df, aes(x = logFC, y = log10_fdr, label = gene, color = sig),
                    inherit.aes = FALSE, size = 2.5, fontface = "bold",
                    max.overlaps = 30, box.padding = 0.4, min.segment.length = 0.2, show.legend = FALSE) +
    scale_color_manual(values = col_vals,
                       labels = c(up = paste0("Up in HB3/4 (n=", n_up, ")"),
                                  down = paste0("Down in HB3/4 (n=", n_down, ")"),
                                  ns = paste0("n.s. (n=", n_ns, ")")), name = NULL) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~FC~(HB3/4~vs~HB1/2/5)), y = expression(-log[10](FDR))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 7.5, color = "grey50"),
          panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 9))
}

deg_results <- mapply(run_deg, files, names(files), SIMPLIFY = FALSE)
plots       <- mapply(make_gene_volcano, deg_results, names(files), SIMPLIFY = FALSE)

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "Gene Volcano — HB3+HB4 vs HB1+HB2+HB5",
                  subtitle = "limma on log2(CPM+1) | positive logFC = higher in HB3/4",
                  theme = theme(plot.title    = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 10, color = "grey40")))

ggsave("gene_volcano_HB34_vs_HB125.pdf", combined, width = 14, height = 11)
ggsave("gene_volcano_HB34_vs_HB125.png", combined, width = 14, height = 11, dpi = 300)
message("Saved!")

for (nm in names(deg_results)) {
  cat("\n===", nm, "===\n")
  res <- deg_results[[nm]]
  res <- res[res$adj.P.Val < FDR_THRESH & abs(res$logFC) > FC_THRESH, ]
  res <- res[order(res$adj.P.Val), ][1:min(15, nrow(res)), c("gene", "logFC", "AveExpr", "adj.P.Val")]
  print(res, row.names = FALSE)
}


######PATHWAY_ANALYSIS_ONLY_USING_SIGNIFICANT_GENES##########
# ============================================================
# ORA Pathway Analysis — significant DEGs (CPM input)
# Fisher's exact test, BH correction, Hallmark gene sets
# ============================================================

library(ggplot2)
library(patchwork)
library(limma)
library(msigdbr)

files <- list(
  "LT-HSC\nLentiNeg"    = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSC\nLentiPos"    = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "Monocytes\nLentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt",
  "Monocytes\nLentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt"
)

FC_THRESH  <- 0.5
FDR_THRESH <- 0.05
N_PATHWAYS <- 15

message("Loading MSigDB Hallmark gene sets (mouse)...")
msig     <- msigdbr(species = "Mus musculus", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

clean_name <- function(x) {
  x <- sub("^HALLMARK_", "", x); x <- gsub("_", " ", x); tools::toTitleCase(tolower(x))
}

run_deg <- function(path, dataset_name) {
  message("\nProcessing: ", dataset_name)
  
  df          <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  hb_cols     <- grep("HB", colnames(df), value = TRUE)
  cpm_mat     <- df[, hb_cols, drop = FALSE]
  
  sample_nums <- sub(".*HB(\\d)$", "\\1", hb_cols)
  group       <- factor(ifelse(sample_nums %in% c("3", "4"), "outlier", "normal"),
                        levels = c("normal", "outlier"))
  
  grp_means <- sapply(levels(group), function(g)
    rowMeans(cpm_mat[, group == g, drop = FALSE]))
  keep    <- rowSums(grp_means >= 1) >= 1
  cpm_mat <- cpm_mat[keep, ]
  
  log_cpm <- log2(cpm_mat + 1)
  design  <- model.matrix(~ group)
  fit     <- eBayes(lmFit(log_cpm, design))
  res     <- topTable(fit, coef = "groupoutlier", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  
  background <- res$gene
  sig_up     <- res$gene[res$adj.P.Val < FDR_THRESH & res$logFC >  FC_THRESH]
  sig_down   <- res$gene[res$adj.P.Val < FDR_THRESH & res$logFC < -FC_THRESH]
  
  message("  Up: ", length(sig_up), "  Down: ", length(sig_down), "  Background: ", length(background))
  list(background = background, up = sig_up, down = sig_down)
}

run_ora <- function(sig_genes, background, pathway_list, min_size = 5) {
  N   <- length(background)
  K   <- length(sig_genes)
  results <- lapply(names(pathway_list), function(pw_name) {
    pw_genes <- intersect(pathway_list[[pw_name]], background)
    M  <- length(pw_genes)
    if (M < min_size) return(NULL)
    k  <- length(intersect(sig_genes, pw_genes))
    p  <- phyper(k - 1, M, N - M, K, lower.tail = FALSE)
    data.frame(pathway = pw_name, p_value = p, n_overlap = k,
               n_sig = K, n_pathway = M, n_bg = N,
               genes_hit = paste(intersect(sig_genes, pw_genes), collapse = ", "),
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, Filter(Negate(is.null), results))
  if (nrow(res) == 0) return(res)
  res$padj          <- p.adjust(res$p_value, method = "BH")
  res$pathway_clean <- clean_name(res$pathway)
  res$fold_enrichment <- (res$n_overlap / res$n_sig) / (res$n_pathway / res$n_bg)
  res[order(res$padj), ]
}

make_ora_plot <- function(ora_up, ora_down, title, n_show = N_PATHWAYS) {
  prep <- function(df, direction) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[df$n_overlap > 0, ]
    if (nrow(df) == 0) return(NULL)
    df <- head(df, n_show); df$direction <- direction; df
  }
  plot_df <- rbind(prep(ora_up, "Up in HB3/4"), prep(ora_down, "Down in HB3/4"))
  if (is.null(plot_df) || nrow(plot_df) == 0)
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No enriched pathways", size=5, color="grey50") + theme_void() + ggtitle(title))
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean,
                                  levels = rev(unique(plot_df$pathway_clean[order(plot_df$direction, plot_df$padj)])))
  plot_df$log10_padj <- -log10(plot_df$padj + 1e-300)
  plot_df$sig_label  <- ifelse(plot_df$padj < 0.05, "*", "")
  
  ggplot(plot_df, aes(x = fold_enrichment, y = pathway_clean, size = n_overlap,
                      color = direction, alpha = log10_padj)) +
    geom_point() +
    geom_text(aes(label = sig_label), size = 5, hjust = -0.8, vjust = 0.8, show.legend = FALSE) +
    facet_wrap(~ direction, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Up in HB3/4" = "#c0392b", "Down in HB3/4" = "#2980b9"), guide = "none") +
    scale_size_continuous(name = "Genes\noverlapping", range = c(2, 8)) +
    scale_alpha_continuous(name = expression(-log[10](padj)), range = c(0.3, 1)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, x = "Fold enrichment", y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "#f0f0f0", color = "#ccc"),
          panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8), legend.position = "right")
}

deg_list <- mapply(run_deg, files, names(files), SIMPLIFY = FALSE)
ora_list <- lapply(deg_list, function(x) list(up = run_ora(x$up, x$background, pathways),
                                              down = run_ora(x$down, x$background, pathways)))
plots    <- mapply(function(ora, title) make_ora_plot(ora$up, ora$down, title),
                   ora_list, names(files), SIMPLIFY = FALSE)

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "ORA — Significant DEGs in HB3+HB4 vs HB1+HB2+HB5",
                  subtitle = sprintf("Fisher's exact test, BH-corrected | FDR < %.2f, |logFC| > %.1f | Hallmark gene sets | limma on log2(CPM+1)",
                                     FDR_THRESH, FC_THRESH),
                  theme = theme(plot.title    = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 10, color = "grey40")))

ggsave("ora_pathway_HB34_vs_HB125.pdf", combined, width = 16, height = 14)
ggsave("ora_pathway_HB34_vs_HB125.png", combined, width = 16, height = 14, dpi = 300)
message("Saved!")

for (nm in names(ora_list)) {
  cat("\n===", nm, "=== UP in HB3/4 ===\n")
  up <- ora_list[[nm]]$up
  if (!is.null(up) && nrow(up) > 0)
    print(head(up[, c("pathway_clean","n_overlap","fold_enrichment","padj")], 10), row.names = FALSE)
  else cat("  No significant pathways\n")
  
  cat("\n===", nm, "=== DOWN in HB3/4 ===\n")
  dn <- ora_list[[nm]]$down
  if (!is.null(dn) && nrow(dn) > 0)
    print(head(dn[, c("pathway_clean","n_overlap","fold_enrichment","padj")], 10), row.names = FALSE)
  else cat("  No significant pathways\n")
}


###HEATMAP_WITH_STRESS_CATEGORY_VOLCANOES###
# ============================================================
# Gene Volcano — stress category colours (CPM input)
# Heat & Stress Response | RNA Handling Artifacts
# limma on log2(CPM+1)
# ============================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
library(limma)
library(msigdbr)

files <- list(
  "LT-HSCs\nLentiNeg"             = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSCs\nLentiPos"             = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "Classical Monocytes\nLentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt",
  "Classical Monocytes\nLentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt"
)

FC_THRESH  <- 0.5
FDR_THRESH <- 0.05
N_LABEL    <- 12

msig_hallmark <- msigdbr(species = "Mus musculus", category = "H")
msig_go_bp    <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")

get_genes <- function(db, set_name) unique(db$gene_symbol[db$gs_name == set_name])

UPR     <- get_genes(msig_hallmark, "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
ROS     <- get_genes(msig_hallmark, "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
HS      <- get_genes(msig_hallmark, "HALLMARK_HEAT_SHOCK")
GO_HEAT <- get_genes(msig_go_bp,   "GOBP_RESPONSE_TO_HEAT")
GO_OX   <- get_genes(msig_go_bp,   "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS")
ALL_HSP <- unique(msig_go_bp$gene_symbol[grepl("^Hsp", msig_go_bp$gene_symbol)])
HEAT_STRESS <- unique(c(HS, UPR, ROS, GO_HEAT, GO_OX, ALL_HSP))

IEG_DISSOCIATION <- c(
  "Fos","Fosb","Fosl1","Fosl2","Jun","Junb","Jund",
  "Egr1","Egr2","Egr3","Egr4",
  "Nr4a1","Nr4a2","Nr4a3",
  "Atf3","Klf4","Klf6",
  "Dusp1","Dusp5","Dusp6",
  "Zfp36","Zfp36l1","Zfp36l2",
  "Ier2","Ier3","Ier5",
  "Rhob","Cyr61","Ctgf","Srf","Myc"
)
COLD_SHOCK   <- c("Cirbp", "Rbm3", "Rbm11")
GO_COLD      <- get_genes(msig_go_bp, "GOBP_RESPONSE_TO_COLD")
RNA_HANDLING <- setdiff(unique(c(IEG_DISSOCIATION, COLD_SHOCK, GO_COLD)), HEAT_STRESS)

CAT_COLS <- c("Heat & Stress" = "#e67e22", "RNA Handling" = "#27ae60",
              "Sig. other"   = "#c0392b", "Sig. down"    = "#2980b9", "n.s." = "#cccccc")

run_deg <- function(path, dataset_name) {
  message("\nProcessing: ", dataset_name)
  df          <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  hb_cols     <- grep("HB", colnames(df), value = TRUE)
  cpm_mat     <- df[, hb_cols, drop = FALSE]
  sample_nums <- sub(".*HB(\\d)$", "\\1", hb_cols)
  group       <- factor(ifelse(sample_nums %in% c("3","4"), "outlier", "normal"), levels = c("normal","outlier"))
  grp_means   <- sapply(levels(group), function(g) rowMeans(cpm_mat[, group == g, drop = FALSE]))
  cpm_mat     <- cpm_mat[rowSums(grp_means >= 1) >= 1, ]
  log_cpm     <- log2(cpm_mat + 1)
  design      <- model.matrix(~ group)
  fit         <- eBayes(lmFit(log_cpm, design))
  res         <- topTable(fit, coef = "groupoutlier", number = Inf, sort.by = "none")
  res$gene    <- rownames(res)
  res
}

annotate_category <- function(df, fc_thresh, fdr_thresh) {
  sig_up  <- df$adj.P.Val < fdr_thresh & df$logFC >  fc_thresh
  sig_down<- df$adj.P.Val < fdr_thresh & df$logFC < -fc_thresh
  in_hs   <- df$gene %in% HEAT_STRESS
  in_rh   <- df$gene %in% RNA_HANDLING
  df$category <- ifelse(!sig_up & !sig_down, "n.s.",
                        ifelse(sig_down,           "Sig. down",
                               ifelse(sig_up & in_hs,     "Heat & Stress",
                                      ifelse(sig_up & in_rh,     "RNA Handling", "Sig. other"))))
  df
}

make_gene_volcano <- function(deg_res, title,
                              fc_thresh  = FC_THRESH, fdr_thresh = FDR_THRESH, n_label = N_LABEL) {
  df           <- annotate_category(deg_res, fc_thresh, fdr_thresh)
  df$log10_fdr <- -log10(df$adj.P.Val + 1e-300)
  n_hs    <- sum(df$category == "Heat & Stress")
  n_rh    <- sum(df$category == "RNA Handling")
  n_other <- sum(df$category == "Sig. other")
  n_down  <- sum(df$category == "Sig. down")
  subtitle <- sprintf("Up — Heat & Stress: %d  |  RNA Handling: %d  |  Other: %d    Down: %d    FDR<%.2f |logFC|>%.1f",
                      n_hs, n_rh, n_other, n_down, fdr_thresh, fc_thresh)
  
  make_labels <- function(cat, n) {
    sub <- df[df$category == cat, ]
    sub[order(sub$adj.P.Val), ][seq_len(min(n, nrow(sub))), ]
  }
  label_df <- rbind(make_labels("Heat & Stress", n_label), make_labels("RNA Handling", ceiling(n_label/2)),
                    make_labels("Sig. other", 6), make_labels("Sig. down", 6))
  if (nrow(label_df) == 0) label_df <- df[order(df$adj.P.Val), ][seq_len(min(n_label, nrow(df))), ]
  
  plot_order <- c("n.s.","Sig. down","Sig. other","RNA Handling","Heat & Stress")
  df$category <- factor(df$category, levels = plot_order)
  df <- df[order(df$category), ]
  
  point_sizes  <- c("n.s." = 0.8, "Sig. down" = 1.6, "Sig. other" = 1.6, "RNA Handling" = 2.2, "Heat & Stress" = 2.5)
  point_alphas <- c("n.s." = 0.3, "Sig. down" = 0.6, "Sig. other" = 0.7, "RNA Handling" = 0.9, "Heat & Stress" = 1.0)
  
  ggplot(df, aes(x = logFC, y = log10_fdr, colour = category, size = category, alpha = category)) +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", colour = "#aaa", linewidth = 0.4) +
    geom_vline(xintercept =  fc_thresh, linetype = "dashed", colour = "#aaa", linewidth = 0.4) +
    geom_vline(xintercept = -fc_thresh, linetype = "dashed", colour = "#aaa", linewidth = 0.4) +
    geom_point() +
    geom_text_repel(data = label_df, aes(x = logFC, y = log10_fdr, label = gene, colour = category),
                    inherit.aes = FALSE, size = 4, fontface = "bold",
                    max.overlaps = 40, box.padding = 0.45, min.segment.length = 0.2, show.legend = FALSE) +
    scale_colour_manual(values = CAT_COLS,
                        labels = c("Heat & Stress" = paste0("Heat & Stress Response (n=", n_hs, ")"),
                                   "RNA Handling"  = paste0("RNA Handling (n=", n_rh, ")"),
                                   "Sig. other"    = paste0("Sig. up, other (n=", n_other, ")"),
                                   "Sig. down"     = paste0("Sig. down (n=", n_down, ")"),
                                   "n.s."          = "n.s."),
                        name = NULL, breaks = c("Heat & Stress","RNA Handling","Sig. other","Sig. down","n.s.")) +
    scale_size_manual(values = point_sizes, guide = "none") +
    scale_alpha_manual(values = point_alphas, guide = "none") +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~FC~~(HB3/4~vs~HB1/2/5)), y = expression(-log[10](FDR))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey50"),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = 12), legend.key.size = unit(0.55, "cm")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_cartesian(ylim = c(0, 4))
}

deg_results <- mapply(run_deg, files, names(files), SIMPLIFY = FALSE)
plots       <- mapply(make_gene_volcano, deg_results, names(files), SIMPLIFY = FALSE)

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title    = "Differential Gene Expression Analysis — HB3+HB4 vs HB1+HB2+HB5",
    subtitle = paste0("Gene sets: MSigDB Hallmark + GO:BP (msigdbr) + all Hsp* family  ·  ",
                      "Dissociation IEGs: van den Brink 2017 (Nat Methods) & Adam 2017 (Cell Rep)  ·  limma on log2(CPM+1)"),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, colour = "grey40")))

ggsave("gene_volcano_HB34_vs_HB125.pdf", combined, width = 15, height = 12)
ggsave("gene_volcano_HB34_vs_HB125.png", combined, width = 15, height = 12, dpi = 300)
message("Saved!")


####better_heatmap####
# ============================================================
# Heatmaps — upregulated heat/stress DEGs (CPM input)
# limma on log2(CPM+1), z-scored per gene for display
# ============================================================

library(limma); library(pheatmap); library(grid); library(msigdbr)

files <- c(
  "LT-HSCs · LentiPos"             = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSCs · LentiNeg"             = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "Classical Monocytes · LentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt",
  "Classical Monocytes · LentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt"
)

FDR_THRESH <- 0.05
FC_THRESH  <- 0.0
MAX_GENES  <- 150

get_heat_stress_genes <- function() {
  get_gs <- function(cat, subcat = NULL, name) {
    args <- list(species = "Mus musculus", category = cat)
    if (!is.null(subcat)) args$subcategory <- subcat
    df <- do.call(msigdbr, args)
    unique(df$gene_symbol[df$gs_name == name])
  }
  genes <- unique(c(get_gs("H", name = "HALLMARK_HEAT_SHOCK"),
                    get_gs("H", name = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
                    get_gs("H", name = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
                    get_gs("C5","GO:BP","GOBP_RESPONSE_TO_HEAT"),
                    get_gs("C5","GO:BP","GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS")))
  go_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
  unique(c(genes, unique(go_bp$gene_symbol[grepl("^Hsp", go_bp$gene_symbol)])))
}

heat_stress_genes <- get_heat_stress_genes()

get_panel_data <- function(path, nm) {
  message("Processing: ", nm)
  df         <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  count_cols <- grep("HB|LB", colnames(df), value = TRUE)
  cpm_all    <- df[, count_cols, drop = FALSE]
  hb_cols    <- grep("HB", count_cols, value = TRUE)
  cpm_hb     <- df[, hb_cols, drop = FALSE]
  
  group <- factor(ifelse(sub(".*HB(\\d)$", "\\1", hb_cols) %in% c("3","4"), "outlier", "normal"),
                  levels = c("normal","outlier"))
  
  # Filter and DEG on HB samples only
  grp_means <- sapply(levels(group), function(g) rowMeans(cpm_hb[, group == g, drop = FALSE]))
  keep      <- rowSums(grp_means >= 1) >= 1
  log_cpm   <- log2(cpm_hb[keep, ] + 1)
  design    <- model.matrix(~ group)
  fit       <- eBayes(lmFit(log_cpm, design))
  res       <- topTable(fit, coef = "groupoutlier", number = Inf, sort.by = "p")
  
  sigs <- rownames(res)[res$adj.P.Val < FDR_THRESH & res$logFC > FC_THRESH &
                          rownames(res) %in% heat_stress_genes]
  if (length(sigs) > MAX_GENES) sigs <- sigs[seq_len(MAX_GENES)]
  message("  Heat/stress upregulated genes: ", length(sigs))
  
  # Display: z-score log2(CPM+1) across ALL samples
  genes <- intersect(sigs, rownames(cpm_all))
  mat_z <- t(scale(t(log2(as.matrix(cpm_all[genes, , drop = FALSE]) + 1))))
  mat_z[is.nan(mat_z)] <- 0
  
  col_order <- c(grep("HB[125]$", colnames(mat_z), value = TRUE),
                 grep("HB[34]$",  colnames(mat_z), value = TRUE),
                 grep("LB",       colnames(mat_z), value = TRUE))
  col_order <- col_order[col_order %in% colnames(mat_z)]
  mat_z     <- mat_z[, col_order, drop = FALSE]
  colnames(mat_z) <- sub("Lenti(Pos|Neg)_", "", colnames(mat_z))
  
  list(mat_z = mat_z, col_order = col_order, n_sig = length(sigs))
}

panels    <- mapply(get_panel_data, files, names(files), SIMPLIFY = FALSE)
pad       <- 0.02
vp_coords <- list(c(0+pad,0.5+pad,0.5-pad,1-pad), c(0.5+pad,0.5+pad,1-pad,1-pad),
                  c(0+pad,0+pad,  0.5-pad,0.5-pad),c(0.5+pad,0+pad,  1-pad,0.5-pad))

draw_panels <- function() {
  grid.newpage()
  for (i in seq_along(panels)) {
    p  <- panels[[i]]; if (p$n_sig == 0) next
    nm <- names(panels)[i]; vp <- vp_coords[[i]]
    gaps <- c(length(grep("HB[125]$", p$col_order)), length(grep("HB[125]$|HB[34]$", p$col_order)))
    pushViewport(viewport(x=vp[1],y=vp[2],width=vp[3]-vp[1],height=vp[4]-vp[2],just=c("left","bottom")))
    print(pheatmap(p$mat_z, color = colorRampPalette(c("#2980b9","white","#c0392b"))(100),
                   cluster_rows = TRUE, cluster_cols = FALSE,
                   show_rownames = TRUE, show_colnames = TRUE,
                   fontsize_row = 20, fontsize_col = 20, fontsize = 24,
                   main = nm, border_color = NA, gaps_col = gaps, silent = TRUE),
          vp = current.viewport())
    popViewport()
  }
}

pdf("heatmaps_heat_stress_up.pdf", width=22, height=20); draw_panels(); dev.off()
png("heatmaps_heat_stress_up.png", width=22, height=20, units="in", res=200); draw_panels(); dev.off()
message("Saved: heatmaps_heat_stress_up.pdf / .png")


###heat_and_Stress_response_score###
# ============================================================
# Heat & Stress Score — all significant genes + Hsp* only
# Score = mean max-scaled log2(CPM+1) across gene set
# Input is CPM: no TMM or cpm() call needed
# ============================================================

library(limma); library(msigdbr)
library(ggplot2); library(dplyr)

files <- c(
  "LT-HSCs · LentiPos"             = "LT_HSCs_DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
  "LT-HSCs · LentiNeg"             = "LT_HSCs_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0.txt",
  "Classical Monocytes · LentiPos" = "Classical_monocytes_DEG_LentiPosHB_vs_LentiPosLB.txt",
  "Classical Monocytes · LentiNeg" = "Classcial_monocytes_DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0 - Copy.txt"
)

FDR_THRESH <- 0.05

get_heat_stress_genes <- function() {
  get_gs <- function(cat, subcat = NULL, name) {
    args <- list(species = "Mus musculus", category = cat)
    if (!is.null(subcat)) args$subcategory <- subcat
    df <- do.call(msigdbr, args)
    unique(df$gene_symbol[df$gs_name == name])
  }
  genes <- unique(c(get_gs("H", name = "HALLMARK_HEAT_SHOCK"),
                    get_gs("H", name = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
                    get_gs("H", name = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
                    get_gs("C5","GO:BP","GOBP_RESPONSE_TO_HEAT"),
                    get_gs("C5","GO:BP","GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS")))
  go_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
  unique(c(genes, grep("^Hsp", unique(go_bp$gene_symbol), value = TRUE)))
}

heat_stress_genes <- get_heat_stress_genes()

run_deg <- function(path, nm) {
  message("DEG: ", nm)
  df          <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  hb_cols     <- grep("HB", colnames(df), value = TRUE)
  cpm_hb      <- df[, hb_cols, drop = FALSE]
  group       <- factor(ifelse(sub(".*HB(\\d)$", "\\1", hb_cols) %in% c("3","4"), "outlier", "normal"),
                        levels = c("normal","outlier"))
  grp_means   <- sapply(levels(group), function(g) rowMeans(cpm_hb[, group == g, drop = FALSE]))
  log_cpm     <- log2(cpm_hb[rowSums(grp_means >= 1) >= 1, ] + 1)
  design      <- model.matrix(~ group)
  fit         <- eBayes(lmFit(log_cpm, design))
  res         <- topTable(fit, coef = "groupoutlier", number = Inf, sort.by = "p")
  rownames(res)[res$adj.P.Val < FDR_THRESH & res$logFC > 0 & rownames(res) %in% heat_stress_genes]
}

sig_per_dataset <- mapply(run_deg, files, names(files), SIMPLIFY = FALSE)
catalogue_all   <- sort(unique(unlist(sig_per_dataset)))
catalogue_hsp   <- sort(grep("^Hsp", catalogue_all, value = TRUE))

message("\nAll heat/stress catalogue: ", length(catalogue_all), " genes")
message("Hsp* catalogue:            ", length(catalogue_hsp), " genes")

# Score: input is already CPM — apply log2(CPM+1) directly, then max-scale per gene
compute_scores <- function(path, nm, gene_set) {
  df         <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  count_cols <- grep("HB|LB", colnames(df), value = TRUE)
  cpm_mat    <- df[, count_cols, drop = FALSE]
  genes      <- intersect(gene_set, rownames(cpm_mat))
  if (length(genes) == 0) return(NULL)
  
  # log2(CPM+1) then max-scale per gene so each gene contributes equally
  # 0 = zero expression, 1 = highest-expressing sample; avoids stretching differences
  mat   <- log2(as.matrix(cpm_mat[genes, , drop = FALSE]) + 1)
  mat01 <- t(apply(mat, 1, function(x) {
    m <- max(x); if (m == 0) return(rep(0, length(x))); x / m
  }))
  score <- colMeans(mat01)
  
  data.frame(
    sample  = sub("Lenti(Pos|Neg)_", "", names(score)),
    score   = score,
    dataset = nm,
    group   = factor(ifelse(grepl("HB[34]$", names(score)), "HB3/4 (outlier)",
                            ifelse(grepl("HB",       names(score)), "HB other", "LB")),
                     levels = c("HB3/4 (outlier)", "HB other", "LB")),
    stringsAsFactors = FALSE
  )
}

scores_all <- bind_rows(mapply(compute_scores, files, names(files),
                               MoreArgs = list(gene_set = catalogue_all), SIMPLIFY = FALSE))
scores_hsp <- bind_rows(mapply(compute_scores, files, names(files),
                               MoreArgs = list(gene_set = catalogue_hsp), SIMPLIFY = FALSE))

scores_all$dataset <- factor(scores_all$dataset, levels = names(files))
scores_hsp$dataset <- factor(scores_hsp$dataset, levels = names(files))

pal <- c("HB3/4 (outlier)" = "#c0392b", "HB other" = "#e07b54", "LB" = "#5b8db8")

make_plot <- function(scores, title) {
  ggplot(scores, aes(x = sample, y = score, fill = group)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    facet_wrap(~ dataset, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(title = title, x = NULL, y = "Mean max-scaled log2(CPM+1)") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          strip.text  = element_text(size = 16, face = "bold"),
          plot.title  = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 14), legend.position = "top",
          panel.grid.major.x = element_blank())
}

p_all <- make_plot(scores_all, sprintf("Heat & Stress Score — all genes (%d genes, ≥1 comparison)", length(catalogue_all)))
p_hsp <- make_plot(scores_hsp, sprintf("Heat & Stress Score — Hsp* only (%d genes)", length(catalogue_hsp)))

ggsave("heat_stress_score_all.pdf", p_all, width = 14, height = 10)
ggsave("heat_stress_score_all.png", p_all, width = 14, height = 10, dpi = 200)
ggsave("heat_stress_score_hsp.pdf", p_hsp, width = 14, height = 10)
ggsave("heat_stress_score_hsp.png", p_hsp, width = 14, height = 10, dpi = 200)
message("Saved: heat_stress_score_all.pdf/png and heat_stress_score_hsp.pdf/png")
message("\nGene catalogue:\n", paste(catalogue_all, collapse = "\n"))


###HB2_HB5_analysis###
setwd("C:/Users/fc809/Downloads/LT-HSCs (1)/LentiPosHB_vs_LentiPosLB")

library(tidyverse)
library(patchwork)
library(ggrepel)

df <- read.delim("DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt",
                 row.names = 1, check.names = FALSE)

sample_cols <- c("LentiPos_HB1","LentiPos_HB2","LentiPos_HB3","LentiPos_HB4","LentiPos_HB5",
                 "LentiPos_LB1","LentiPos_LB2","LentiPos_LB3","LentiPos_LB4")

# Input is CPM — use directly
expr  <- df[, sample_cols]
short <- gsub("LentiPos_", "", sample_cols)
colnames(expr) <- short

group    <- ifelse(grepl("HB", short), "HB", "LB")
hb_color <- "#E74C3C"; lb_color <- "#3498DB"
col_scale <- c(HB = hb_color, LB = lb_color)

panel_theme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold",size=12,margin=margin(b=2)),
        plot.subtitle = element_text(size=10,color="grey50",margin=margin(b=6)),
        plot.margin = margin(t=6,r=10,b=6,l=6), legend.position="top",
        legend.margin = margin(t=0,b=2), legend.key.size = unit(0.7,"lines"),
        legend.text = element_text(size=10), legend.box.spacing = unit(2,"pt"),
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        panel.grid = element_line(color="grey93"))

# ── Panel A: PCA on log2(CPM+1) ──────────────────────────────
expr_log <- log2(t(expr) + 1)
pca_res  <- prcomp(scale(expr_log), center = FALSE, scale. = FALSE)

pca_df <- data.frame(sample = short, group = group,
                     PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                     outlier = short %in% c("HB2","HB5"))
var_exp <- round(summary(pca_res)$importance[2,1:2]*100, 1)

pA <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(aes(size = outlier, shape = outlier), stroke = 1.2) +
  geom_text_repel(aes(label = sample), size = 3.5, show.legend = FALSE,
                  box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf) +
  scale_color_manual(values = col_scale, name = NULL) +
  scale_size_manual(values  = c(`FALSE`=2.5, `TRUE`=3.5), guide="none") +
  scale_shape_manual(values = c(`FALSE`=16,  `TRUE`=21),  guide="none") +
  labs(title = "Sample PCA", subtitle = "Pseudo-bulk expression log2(CPM+1)",
       x = paste0("PC1 (", var_exp[1], "% variance)"),
       y = paste0("PC2 (", var_exp[2], "% variance)")) +
  panel_theme

# ── Panel B: Gene dropout ─────────────────────────────────────
zero_df <- data.frame(sample = short, group = group,
                      n_zeros = colSums(expr == 0), outlier = short %in% c("HB2","HB5"))

pB <- ggplot(zero_df, aes(x = fct_inorder(sample), y = n_zeros, fill = group)) +
  geom_col(aes(color = outlier), linewidth = 0.9) +
  geom_text(aes(label = ifelse(n_zeros > 0, n_zeros, "")), vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = col_scale, name = NULL) +
  scale_color_manual(values = c(`FALSE`=NA,`TRUE`="black"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Gene dropout per sample", subtitle = "Number of genes with zero expression", x = NULL, y = "Count") +
  panel_theme +
  theme(axis.text.x = element_text(angle=45,hjust=1), panel.grid.major.x = element_blank(), legend.position="none")

# ── Panel C: CPM fraction by gene set ────────────────────────
# Since input is CPM, use it directly as proportional measure
prolif_genes <- c("Mki67","Pclaf","Cdc20","Birc5","Ube2c","Cdca8","Cks1b","Cdc6")
oxphos_genes <- c("Ndufv1","Sdha","Cox4i1","Atp5b","Cycs","Cox5a","Uqcrc2","Ndufs1")

found_prolif <- intersect(prolif_genes, rownames(expr))
found_oxphos <- intersect(oxphos_genes, rownames(expr))
total_cpm    <- colSums(expr)

freq_df <- data.frame(
  sample         = short,
  group          = group,
  Proliferation  = colSums(expr[found_prolif, ]) / total_cpm * 100,
  `Mito. OXPHOS` = colSums(expr[found_oxphos, ]) / total_cpm * 100,
  check.names    = FALSE
) %>% pivot_longer(cols = c("Proliferation","Mito. OXPHOS"), names_to = "gene_set", values_to = "pct")

pC <- ggplot(freq_df, aes(x = fct_inorder(sample), y = pct, fill = gene_set)) +
  geom_col(aes(color = sample == "HB5"), position = position_dodge(width = 0.75), width = 0.7, linewidth = 0.9) +
  annotate("text", x = "HB5", y = max(freq_df$pct)*1.1,
           label = "elevated transcript\nburden", size = 3, fontface = "italic", hjust = 0.5, color = "grey35") +
  scale_fill_manual(values = c("Proliferation" = hb_color, "Mito. OXPHOS" = lb_color), name = NULL) +
  scale_color_manual(values = c(`FALSE`=NA,`TRUE`="black"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "HB5 biological signature",
       subtitle = "Transcript frequency (% of total CPM) by gene set", x = NULL, y = "% of total CPM") +
  panel_theme +
  theme(axis.text.x = element_text(angle=45,hjust=1), panel.grid.major.x = element_blank())

combined <- pA + pB + pC +
  plot_annotation(title = "HB2 and HB5 outlier analysis — LentiPos pseudo-bulk RNA-seq",
                  theme = theme(plot.title = element_text(face="bold",size=13),
                                plot.margin = margin(t=8,r=8,b=4,l=8)))

ggsave("outlier_analysis_combined.pdf", combined, width = 16, height = 6.5)
ggsave("outlier_analysis_combined.png", combined, width = 16, height = 6.5, dpi = 160)
message("Saved: outlier_analysis_combined.pdf / .png")

