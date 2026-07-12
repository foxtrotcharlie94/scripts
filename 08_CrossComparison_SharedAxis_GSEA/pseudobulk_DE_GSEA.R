# ============================================================
# Pseudobulk DE + GSEA — fully self-contained from 3 raw inputs:
#   1. count_matrix.txt
#        TSV: genes x cells, first column = gene names,
#        header row = cell barcodes (format "HB1_AAACC...")
#   2. sample_and_clusterNumber_for_each_cell_res0.5.txt
#        columns: cells, sample, harmony_clusters_res0.5
#   3. barcode_counts_raw_with_sort.csv  (optional)
#        per-cell sort metadata
#
# Pool one or more clusters into `target_label` and run 4 comparisons:
#   1. HB:       LentiPos vs LentiNeg   (paired by donor)
#   2. LB:       LentiPos vs LentiNeg   (paired by donor)
#   3. LentiPos: HB vs LB               (unpaired, batch covariate)
#   4. LentiNeg: HB vs LB               (unpaired, batch covariate)
#
# Per comparison: PCA panel, BCV, MD, volcano, 3 sorted DE CSVs,
# GSEA Hallmark + GO_BP each ranked by logFC and logFC*-log10(P).
# Group-level summary at the end (2 bar plots).
# ============================================================

suppressPackageStartupMessages({
  library(edgeR); library(data.table); library(Matrix)
  library(dplyr); library(ggplot2); library(tidyr)
  library(fgsea); library(msigdbr); library(ggrepel)
  library(patchwork); library(RColorBrewer)
})
# Resource budget: use most of the machine.
n_cores <- max(1, parallel::detectCores() - 1)
data.table::setDTthreads(n_cores)
RhpcBLASctl::blas_set_num_threads(n_cores)
RhpcBLASctl::omp_set_num_threads(n_cores)
options(future.globals.maxSize = 64 * 1024^3)  # 64 GB ceiling for future-backed steps
cat(sprintf("Cores detected: %d | using: %d | data.table threads: %d\n",
            parallel::detectCores(), n_cores, getDTthreads()))

# ══════════════════════════════════════════════════════════════════════════════
# Settings — set target_clusters and target_label per scenario
# ══════════════════════════════════════════════════════════════════════════════
in_mat_path  <- "C:/Users/fc809/Downloads/count_matrix.txt"
in_meta_path <- "C:/Users/fc809/Downloads/sample_and_clusterNumber_for_each_cell_res0.5.txt"
in_sort_csv  <- "C:/Users/fc809/Downloads/barcode_counts_raw_with_sort.csv"  # optional

target_clusters <- c("0")                # one or more cluster IDs to pool
target_label    <- "ClassicalMonocytes"   # name used in folders/titles/legends

base_out_dir <- file.path("C:/Users/fc809/Downloads/pseudobulk_DE_GSEA", target_label)
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)

# Lenti classification (ambiguous in-between cells are dropped)
lenti_gene_name <- "LentiAll"
THR_POS <- 5e-5
THR_NEG <- 1e-5

# Filtering & QC parameters
apply_detection_filter <- TRUE
min_pct_per_sample_DE  <- 0.05
min_pct_per_sample_PCA <- 0.25
min_cells_per_pb       <- 10

# GSEA collections + plot rules
gsea_species <- "Mus musculus"
gsea_specs   <- list(
  Hallmark = list(cat = "H",  subcat = NULL,    plot_rule = "sig_only"),
  GO_BP    = list(cat = "C5", subcat = "GO:BP", plot_rule = "top_per_direction",
                  top_per_dir = 20)
)
gsea_minSize     <- 15
gsea_maxSize     <- 500
gsea_nPermSimple <- 10000
sig_padj         <- 0.05

ranking_specs <- list(
  logFC   = list(col = "logFC",           display = "logFC"),
  logFCxP = list(col = "logFC_x_negLogP", display = "logFC * -log10(P)")
)
volcano_top_n <- 30

batch_map <- c("HB1"="B1","HB2"="B1","LB1"="B1","LB2"="B1",
               "HB3"="B2","HB4"="B2","LB3"="B2","LB4"="B2",
               "HB5"="B3")

get_comp_dir <- function(label) {
  d <- file.path(base_out_dir, label)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# ══════════════════════════════════════════════════════════════════════════════
# Load: cluster CSV + selected columns of count_matrix → counts + meta
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading cluster assignments...\n")
meta_cl <- fread(in_meta_path)
setnames(meta_cl, c("cells", "harmony_clusters_res0.5"), c("cell", "cluster"))
meta_cl[, cell := as.character(cell)]
meta_cl[, cluster := as.character(cluster)]
cat(sprintf("  %d cells across %d clusters\n",
            nrow(meta_cl), length(unique(meta_cl$cluster))))

if (!all(target_clusters %in% meta_cl$cluster))
  stop("Target cluster(s) not in metadata: ",
       paste(setdiff(target_clusters, meta_cl$cluster), collapse = ","))

cells_target <- meta_cl[cluster %in% target_clusters, cell]
cat(sprintf("Target clusters (%s): %d cells\n",
            paste(target_clusters, collapse=","), length(cells_target)))

cat("\nReading count_matrix header...\n")
hdr <- names(fread(in_mat_path, nrows = 0))
gene_col <- hdr[1]
cells_in_mat <- hdr[-1]
cells_to_read <- intersect(cells_target, cells_in_mat)
cat(sprintf("  %d / %d target cells present in count_matrix\n",
            length(cells_to_read), length(cells_target)))
if (length(cells_to_read) == 0) stop("No target cells found in count_matrix.")

cat("Reading count_matrix (selected cells only)...\n")
t0 <- Sys.time()
dt <- fread(in_mat_path, select = c(gene_col, cells_to_read))
cat(sprintf("  Done in %.1f min | %d rows x %d cols\n",
            as.numeric(Sys.time() - t0, units = "mins"), nrow(dt), ncol(dt)))

genes <- as.character(dt[[1]])
mat <- as.matrix(dt[, -1, with = FALSE]); rownames(mat) <- genes
counts <- as(mat, "CsparseMatrix"); rm(dt, mat); gc()
cat(sprintf("Counts matrix: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# Lenti classification from the lenti gene row
if (!lenti_gene_name %in% rownames(counts))
  stop("Lenti gene '", lenti_gene_name, "' not found in counts matrix.")
lenti_umis <- as.numeric(counts[lenti_gene_name, ])
total_umis <- as.numeric(Matrix::colSums(counts))
ratio      <- lenti_umis / pmax(total_umis, 1)
genotype   <- ifelse(ratio >= THR_POS, "LentiPos",
              ifelse(ratio <= THR_NEG, "LentiNeg", NA_character_))
cat(sprintf("LentiPos %d | LentiNeg %d | ambiguous (dropped) %d\n",
            sum(genotype == "LentiPos", na.rm = TRUE),
            sum(genotype == "LentiNeg", na.rm = TRUE),
            sum(is.na(genotype))))

# Drop LentiAll from counts (not a host gene) and drop ambiguous cells
counts <- counts[rownames(counts) != lenti_gene_name, , drop = FALSE]
keep   <- !is.na(genotype)
counts <- counts[, keep, drop = FALSE]
cells_kept <- colnames(counts)

# Build metadata
meta <- data.frame(row.names = cells_kept, cell = cells_kept, stringsAsFactors = FALSE)
meta$sample      <- sub("_.*$", "", meta$cell)
meta$condition   <- ifelse(grepl("^HB", meta$sample), "HB", "LB")
meta$batch       <- batch_map[meta$sample]
meta$cluster     <- meta_cl[match(meta$cell, meta_cl$cell), cluster]
meta$lenti_umis  <- lenti_umis[keep]
meta$total_umis  <- total_umis[keep]
meta$lenti_ratio <- ratio[keep]
meta$genotype    <- genotype[keep]

# Optional sort-file join
if (file.exists(in_sort_csv)) {
  sort_dt <- fread(in_sort_csv)
  cc <- intersect(c("cell","cells","barcode","Barcode","cell_id","CellID","Cell","cell_name"),
                  colnames(sort_dt))[1]
  if (!is.na(cc)) {
    sort_dt[[cc]] <- as.character(sort_dt[[cc]])
    sj <- sort_dt[match(cells_kept, sort_dt[[cc]]), ]
    sj[[cc]] <- NULL
    setnames(sj, paste0("sort_", names(sj)))
    meta <- cbind(meta, as.data.frame(sj))
    cat(sprintf("Joined %d sort columns; matched %d / %d cells.\n",
                ncol(sj), sum(!is.na(sj[[1]])), length(cells_kept)))
  }
}

cat(sprintf("\nGroup '%s' (clusters %s): %d cells\n",
            target_label, paste(target_clusters, collapse = ","), nrow(meta)))
cat("Cells per sample x genotype:\n"); print(table(meta$sample, meta$genotype))

# ══════════════════════════════════════════════════════════════════════════════
# Helpers (operate on (counts, meta) — no Seurat)
# ══════════════════════════════════════════════════════════════════════════════

subset_set <- function(counts, meta, idx) {
  list(counts = counts[, idx, drop = FALSE], meta = meta[idx, , drop = FALSE])
}

build_pb <- function(counts, group_vec, features = NULL) {
  if (!is.null(features)) counts <- counts[features, , drop = FALSE]
  groups <- unique(group_vec)
  M <- sparseMatrix(
    i = seq_along(group_vec), j = match(group_vec, groups), x = 1,
    dims = c(length(group_vec), length(groups)),
    dimnames = list(NULL, groups))
  pb <- as.matrix(counts %*% M); storage.mode(pb) <- "integer"; pb
}

detect_pct_per_sample <- function(counts, sample_vec, min_pct) {
  samples <- unique(sample_vec)
  det_mat <- sapply(samples, function(s) {
    cols <- which(sample_vec == s)
    if (length(cols) == 0) return(rep(FALSE, nrow(counts)))
    as.numeric(Matrix::rowSums(counts[, cols, drop = FALSE] > 0)) / length(cols) >= min_pct
  })
  rownames(counts)[rowSums(det_mat) == length(samples)]
}

get_msig <- function(species, cat, subcat) {
  args_new <- if (is.null(subcat)) list(species = species, collection = cat)
              else list(species = species, collection = cat, subcollection = subcat)
  args_old <- if (is.null(subcat)) list(species = species, category = cat)
              else list(species = species, category = cat, subcategory = subcat)
  res <- tryCatch(do.call(msigdbr, args_new), error = function(e) NULL)
  if (is.null(res)) res <- do.call(msigdbr, args_old)
  sym_col <- intersect(c("gene_symbol","ensembl_gene","gene_symbol_or_ensembl"),
                       colnames(res))[1]
  if (is.na(sym_col)) sym_col <- "gene_symbol"
  list(df = res, sym_col = sym_col)
}

parse_pb_meta <- function(col_names, paired) {
  if (paired)
    data.frame(col = col_names, sample = sub("__.*$", "", col_names),
               genotype = sub("^.*__", "", col_names), stringsAsFactors = FALSE)
  else
    data.frame(col = col_names, sample = col_names, stringsAsFactors = FALSE)
}

do_pca <- function(pb, col_meta, color_var, title, transform = "log_cpm") {
  pb <- pb[rowSums(pb) > 0, , drop = FALSE]
  if (nrow(pb) < 10 || ncol(pb) < 3)
    return(ggplot() + theme_void() + ggtitle(paste(title, "\n(too few features/samples)")))
  y <- DGEList(counts = pb); y <- normLibSizes(y)
  log_cpm <- cpm(y, log = TRUE, prior.count = 3)
  if (transform == "rank") { mat <- t(apply(log_cpm, 1, rank)); do_scale <- FALSE
  } else { mat <- log_cpm; do_scale <- TRUE }
  if (nrow(mat) > 2000) { vars <- apply(mat, 1, var); mat <- mat[order(-vars)[1:2000], ] }
  pca <- prcomp(t(mat), scale. = do_scale)
  ve  <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  d <- merge(data.frame(col = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2]),
             col_meta, by = "col", all.x = TRUE)
  d$batch <- factor(batch_map[d$sample])
  ggplot(d, aes(PC1, PC2, color = .data[[color_var]], shape = batch)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text_repel(aes(label = sample), size = 3, color = "grey30",
                    max.overlaps = 30, segment.color = "grey70") +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2]), title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "plain"))
}

make_pca_panel <- function(set, group_col, color_col, paired, label) {
  comp_dir <- get_comp_dir(label)
  group_vec <- set$meta[[group_col]]
  pb1 <- build_pb(set$counts, group_vec)
  cm  <- parse_pb_meta(colnames(pb1), paired)
  cm$condition <- ifelse(grepl("^HB", cm$sample), "HB", "LB")

  g25 <- detect_pct_per_sample(set$counts, set$meta$sample, min_pct_per_sample_PCA)
  pb2 <- if (length(g25) >= 10) build_pb(set$counts, group_vec, features = g25) else NULL

  mean_per_sample <- aggregate(total_umis ~ sample, data = set$meta, FUN = mean)
  x_thr <- 0.9 * min(mean_per_sample$total_umis)
  keep_d <- set$meta$total_umis >= x_thr
  set_d  <- subset_set(set$counts, set$meta, keep_d)
  pb3 <- if (ncol(set_d$counts) >= 50) build_pb(set_d$counts, set_d$meta[[group_col]]) else NULL
  cm3 <- if (!is.null(pb3)) {
    tmp <- parse_pb_meta(colnames(pb3), paired)
    tmp$condition <- ifelse(grepl("^HB", tmp$sample), "HB", "LB"); tmp
  } else NULL

  one_or_blank <- function(pb, cm_use, base_title, transform) {
    if (is.null(pb) || is.null(cm_use))
      return(ggplot() + theme_void() + ggtitle(paste(base_title, "(skipped)")))
    suffix <- if (transform == "rank") " — rank" else " — log-CPM"
    do_pca(pb, cm_use, color_col, paste0(base_title, suffix), transform = transform)
  }

  t1 <- sprintf("All genes (%d)", if (is.null(pb1)) 0 else nrow(pb1))
  t2 <- sprintf("Genes >=%.0f%% per sample (%d)", 100*min_pct_per_sample_PCA, length(g25))
  t3 <- sprintf("Cells with total >= %.0f UMIs (%d cells)",
                x_thr, if (is.null(pb3)) 0 else ncol(set_d$counts))
  p1 <- one_or_blank(pb1, cm,  t1, "log_cpm")
  p2 <- one_or_blank(pb2, cm,  t2, "log_cpm")
  p3 <- one_or_blank(pb3, cm3, t3, "log_cpm")
  p4 <- one_or_blank(pb1, cm,  t1, "rank")
  p5 <- one_or_blank(pb2, cm,  t2, "rank")
  p6 <- one_or_blank(pb3, cm3, t3, "rank")
  panel <- ((p1 | p2 | p3) / (p4 | p5 | p6)) +
    plot_annotation(
      title    = sprintf("PCA panel — %s | %s", target_label, label),
      subtitle = sprintf("Top: log-CPM | Bottom: per-gene rank | clusters: %s",
                         paste(target_clusters, collapse = ",")))
  ggsave(file.path(comp_dir, sprintf("PCA_%s.png", label)), panel,
         width = 18, height = 13, dpi = 150)
}

finalize_DE <- function(qlf, y, label, contrast_desc = "") {
  comp_dir <- get_comp_dir(label)
  tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tt$gene <- rownames(tt)
  pP <- pmax(tt$PValue, .Machine$double.xmin)
  pF <- pmax(tt$FDR,    .Machine$double.xmin)
  tt$signed_negLogP    <- sign(tt$logFC) * -log10(pP)
  tt$signed_negLogFDR  <- sign(tt$logFC) * -log10(pF)
  tt$logFC_x_negLogP   <- tt$logFC * -log10(pP)
  tt$logFC_x_negLogFDR <- tt$logFC * -log10(pF)
  tt <- tt[, c("gene","logFC","logCPM","F","PValue","FDR",
               "signed_negLogP","signed_negLogFDR",
               "logFC_x_negLogP","logFC_x_negLogFDR")]
  fwrite(tt[order(tt$PValue), ],                  file.path(comp_dir, sprintf("DE_%s_byPValue.csv", label)))
  fwrite(tt[order(-abs(tt$logFC_x_negLogFDR)), ], file.path(comp_dir, sprintf("DE_%s_by_logFCxnegLogFDR.csv", label)))
  fwrite(tt[order(-abs(tt$logFC)), ],             file.path(comp_dir, sprintf("DE_%s_by_absLogFC.csv", label)))
  cat(sprintf("Sig @ FDR<0.05: %d  (up %d, down %d)\n",
              sum(tt$FDR < 0.05), sum(tt$FDR < 0.05 & tt$logFC > 0),
              sum(tt$FDR < 0.05 & tt$logFC < 0)))
  png(file.path(comp_dir, sprintf("BCV_%s.png", label)),
      width = 6, height = 5, units = "in", res = 150)
  plotBCV(y, main = sprintf("BCV — %s | %s", target_label, label)); dev.off()
  png(file.path(comp_dir, sprintf("MD_%s.png", label)),
      width = 6, height = 5, units = "in", res = 150)
  plotMD(qlf, status = decideTests(qlf),
         main = sprintf("MD — %s | %s", target_label, label))
  abline(h = 0, col = "grey60"); dev.off()
  tt
}

run_edgeR_paired <- function(set, label) {
  cat(sprintf("\n========== %s ==========\n", label))
  cat("Cells:", nrow(set$meta), "\n")
  set$meta$pb_group <- paste(set$meta$sample, set$meta$genotype, sep = "__")
  make_pca_panel(set, "pb_group", "genotype", paired = TRUE, label = label)

  gene_keep <- if (apply_detection_filter)
    detect_pct_per_sample(set$counts, set$meta$sample, min_pct_per_sample_DE)
    else rownames(set$counts)
  if (apply_detection_filter)
    cat(sprintf("DE detection >= %.0f%%: %d / %d genes\n",
                100*min_pct_per_sample_DE, length(gene_keep), nrow(set$counts)))

  gs <- table(set$meta$pb_group)
  small <- names(gs[gs < min_cells_per_pb])
  if (length(small) > 0) {
    cat("Dropping pb groups <", min_cells_per_pb, "cells:", paste(small, collapse=", "), "\n")
    set <- subset_set(set$counts, set$meta, !set$meta$pb_group %in% small)
  }
  pb <- build_pb(set$counts, set$meta$pb_group, features = gene_keep)
  cm <- parse_pb_meta(colnames(pb), paired = TRUE); rownames(cm) <- cm$col

  d <- table(cm$sample, cm$genotype)
  paired_donors <- rownames(d)[rowSums(d > 0) == 2]
  cat("Paired donors:", length(paired_donors), "/", nrow(d), "\n")
  if (length(paired_donors) < 2) { cat("Skipped (insufficient pairs).\n"); return(NULL) }
  cm_p <- cm[cm$sample %in% paired_donors, , drop = FALSE]
  pb_p <- pb[, rownames(cm_p), drop = FALSE]
  cm_p$sample   <- factor(cm_p$sample)
  cm_p$genotype <- factor(cm_p$genotype, levels = c("LentiNeg", "LentiPos"))

  y <- DGEList(counts = pb_p, samples = cm_p, group = cm_p$genotype)
  keep <- filterByExpr(y, group = cm_p$genotype)
  cat("After filterByExpr:", sum(keep), "/", length(keep), "genes\n")
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
  design <- model.matrix(~ sample + genotype, data = cm_p)
  y   <- estimateDisp(y, design); fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "genotypeLentiPos")
  list(tt = finalize_DE(qlf, y, label), y = y, qlf = qlf,
       contrast_desc = "LentiPos vs LentiNeg (positive logFC = up in LentiPos)")
}

run_edgeR_unpaired <- function(set, label) {
  cat(sprintf("\n========== %s ==========\n", label))
  cat("Cells:", nrow(set$meta), "\n")
  cat("Cells per sample:\n"); print(table(set$meta$sample))
  make_pca_panel(set, "sample", "condition", paired = FALSE, label = label)

  gene_keep <- if (apply_detection_filter)
    detect_pct_per_sample(set$counts, set$meta$sample, min_pct_per_sample_DE)
    else rownames(set$counts)
  if (apply_detection_filter)
    cat(sprintf("DE detection >= %.0f%%: %d / %d genes\n",
                100*min_pct_per_sample_DE, length(gene_keep), nrow(set$counts)))

  ss <- table(set$meta$sample)
  small <- names(ss[ss < min_cells_per_pb])
  if (length(small) > 0) {
    cat("Dropping samples <", min_cells_per_pb, "cells:", paste(small, collapse=", "), "\n")
    set <- subset_set(set$counts, set$meta, !set$meta$sample %in% small)
  }
  pb <- build_pb(set$counts, set$meta$sample, features = gene_keep)
  cm <- parse_pb_meta(colnames(pb), paired = FALSE); rownames(cm) <- cm$col
  cm$condition <- ifelse(grepl("^HB", cm$sample), "HB", "LB")
  cm$batch     <- batch_map[cm$sample]
  cat("Pseudobulk samples:\n"); print(cm)

  if (sum(cm$condition == "HB") < 2 || sum(cm$condition == "LB") < 2) {
    cat("Need >= 2 samples per condition — skipping\n"); return(NULL)
  }
  cm$batch     <- factor(cm$batch)
  cm$condition <- factor(cm$condition, levels = c("LB", "HB"))
  y <- DGEList(counts = pb, samples = cm, group = cm$condition)
  keep <- filterByExpr(y, group = cm$condition)
  cat("After filterByExpr:", sum(keep), "/", length(keep), "genes\n")
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
  design <- if (nlevels(cm$batch) > 1) model.matrix(~ batch + condition, data = cm)
            else                      model.matrix(~ condition, data = cm)
  y   <- estimateDisp(y, design); fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "conditionHB")
  list(tt = finalize_DE(qlf, y, label), y = y, qlf = qlf,
       contrast_desc = "HB vs LB (positive logFC = up in HB)")
}

plot_volcano <- function(tt, label, contrast_desc, top_n = volcano_top_n) {
  pP <- pmax(tt$PValue, .Machine$double.xmin)
  d <- tt; d$neg_log10P <- -log10(pP)
  d$direction <- ifelse(d$FDR < 0.05 & d$logFC >  0.5, "Up",
                 ifelse(d$FDR < 0.05 & d$logFC < -0.5, "Down", "ns"))
  top_lab <- head(d[order(-abs(d$logFC_x_negLogFDR)), ], top_n)
  ggplot(d, aes(logFC, neg_log10P)) +
    geom_point(aes(color = direction), size = 0.9, alpha = 0.6) +
    scale_color_manual(values = c("Up"="firebrick","Down"="steelblue","ns"="grey80")) +
    geom_text_repel(data = top_lab, aes(label = gene),
                    size = 3, max.overlaps = 40, segment.color = "grey70") +
    labs(x = "log2 FC", y = "-log10(P-value)",
         title = sprintf("Volcano — %s | %s", target_label, label),
         subtitle = sprintf("%s | top %d labeled", contrast_desc, top_n),
         color = NULL) +
    theme_bw(base_size = 12)
}

run_one_gsea <- function(tt, collection_name, cat_code, subcat_code, ranking_col) {
  msig <- get_msig(gsea_species, cat_code, subcat_code)
  pathways <- split(msig$df[[msig$sym_col]], msig$df$gs_name)
  ranks <- tt[[ranking_col]]; names(ranks) <- tt$gene
  ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  set.seed(42)
  fres <- fgsea(pathways = pathways, stats = ranks,
                minSize = gsea_minSize, maxSize = gsea_maxSize,
                nPermSimple = gsea_nPermSimple)
  fres[order(padj)]
}

plot_gsea <- function(fres, label, contrast_desc, collection_name, spec) {
  if (spec$plot_rule == "sig_only") {
    d <- fres[!is.na(padj) & padj < sig_padj]
    title_suffix <- sprintf("(all sig at padj<%g, n=%d)", sig_padj, nrow(d))
  } else if (spec$plot_rule == "top_per_direction") {
    top_e <- head(fres[NES > 0][order(padj)], spec$top_per_dir)
    top_d <- head(fres[NES < 0][order(padj)], spec$top_per_dir)
    d <- rbind(top_e, top_d)
    title_suffix <- sprintf("(top %d enriched + top %d depleted)",
                            spec$top_per_dir, spec$top_per_dir)
  } else { d <- fres; title_suffix <- "" }
  if (nrow(d) == 0)
    return(ggplot() + theme_void() +
           ggtitle(sprintf("GSEA %s — %s: nothing to plot", collection_name, label)))
  d <- d[order(-NES)]
  d[, pathway_label := factor(pathway, levels = rev(unique(pathway)))]
  d[, is_sig := !is.na(padj) & padj < sig_padj]
  ggplot(d, aes(NES, pathway_label,
                fill = -log10(pmax(padj, .Machine$double.xmin)),
                color = is_sig)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_viridis_c(name = "-log10(padj)") +
    scale_color_manual(values = c("TRUE"="black","FALSE"=NA), guide = "none") +
    labs(x = "NES", y = NULL,
         title = sprintf("GSEA %s — %s | %s %s",
                         collection_name, target_label, label, title_suffix),
         subtitle = sprintf("%s | Outlined: padj<%g", contrast_desc, sig_padj)) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7))
}

run_all_for <- function(de_res, label) {
  if (is.null(de_res)) return(invisible(NULL))
  comp_dir <- get_comp_dir(label)
  ggsave(file.path(comp_dir, sprintf("volcano_%s.png", label)),
         plot_volcano(de_res$tt, label, de_res$contrast_desc),
         width = 8, height = 7, dpi = 150)
  for (coll_name in names(gsea_specs)) {
    spec <- gsea_specs[[coll_name]]
    for (rank_name in names(ranking_specs)) {
      rank_col  <- ranking_specs[[rank_name]]$col
      rank_disp <- ranking_specs[[rank_name]]$display
      fres <- run_one_gsea(de_res$tt, coll_name, spec$cat, spec$subcat, rank_col)
      fsave <- copy(fres); fsave[, leadingEdge := sapply(leadingEdge, paste, collapse=";")]
      fwrite(fsave, file.path(comp_dir,
             sprintf("GSEA_%s_%s_rankBy_%s.csv", coll_name, label, rank_name)))
      cat(sprintf("GSEA %s %s [rank: %s] — sig @ padj<%g: %d / %d\n",
                  coll_name, label, rank_disp, sig_padj,
                  sum(fres$padj < sig_padj, na.rm = TRUE), nrow(fres)))
      n_rows <- if (spec$plot_rule == "sig_only") sum(fres$padj < sig_padj, na.rm = TRUE)
                else 2 * spec$top_per_dir
      h <- max(6, min(0.25 * n_rows + 2, 30))
      contrast_with_rank <- sprintf("%s | rank: %s", de_res$contrast_desc, rank_disp)
      ggsave(file.path(comp_dir,
             sprintf("GSEA_%s_%s_rankBy_%s.png", coll_name, label, rank_name)),
             plot_gsea(fres, label, contrast_with_rank, coll_name, spec),
             width = 11, height = h, dpi = 150, limitsize = FALSE)
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Run the four comparisons
# ══════════════════════════════════════════════════════════════════════════════
full_set <- list(counts = counts, meta = meta)
rm(counts); gc()

set_HB <- subset_set(full_set$counts, full_set$meta, full_set$meta$condition == "HB")
set_LB <- subset_set(full_set$counts, full_set$meta, full_set$meta$condition == "LB")
de_HB <- run_edgeR_paired(set_HB, "HB_PosVsNeg"); rm(set_HB); gc()
de_LB <- run_edgeR_paired(set_LB, "LB_PosVsNeg"); rm(set_LB); gc()
run_all_for(de_HB, "HB_PosVsNeg")
run_all_for(de_LB, "LB_PosVsNeg")

set_LP <- subset_set(full_set$counts, full_set$meta, full_set$meta$genotype == "LentiPos")
set_LN <- subset_set(full_set$counts, full_set$meta, full_set$meta$genotype == "LentiNeg")
de_LP <- run_edgeR_unpaired(set_LP, "LentiPos_HBvsLB"); rm(set_LP); gc()
de_LN <- run_edgeR_unpaired(set_LN, "LentiNeg_HBvsLB"); rm(set_LN); gc()
run_all_for(de_LP, "LentiPos_HBvsLB")
run_all_for(de_LN, "LentiNeg_HBvsLB")

# ══════════════════════════════════════════════════════════════════════════════
# Per-group summary: two bar plots
# ══════════════════════════════════════════════════════════════════════════════
comp_order <- c("HB_PosVsNeg", "LB_PosVsNeg", "LentiPos_HBvsLB", "LentiNeg_HBvsLB")
de_rows  <- list(); gse_rows <- list()
for (comp in comp_order) {
  cd <- file.path(base_out_dir, comp)
  if (!dir.exists(cd)) next
  de_fp <- file.path(cd, sprintf("DE_%s_byPValue.csv", comp))
  if (file.exists(de_fp)) {
    d <- fread(de_fp)
    de_rows[[length(de_rows) + 1]] <- data.frame(
      comparison = comp, n_genes_FDR025 = sum(d$FDR < 0.25, na.rm = TRUE))
  }
  gsea_files <- list.files(cd, pattern = "^GSEA_.*_rankBy_.*\\.csv$", full.names = TRUE)
  sig_keys <- character()
  for (fp in gsea_files) {
    bn <- basename(fp)
    m  <- regmatches(bn, regexec(sprintf("^GSEA_(.+)_%s_rankBy_(.+)\\.csv$", comp), bn))[[1]]
    coll <- if (length(m) >= 3) m[2] else "unknown"
    d <- fread(fp)
    s <- d[!is.na(padj) & padj < sig_padj]
    if (nrow(s) > 0) sig_keys <- union(sig_keys, paste(coll, s$pathway, sep = "|"))
  }
  gse_rows[[length(gse_rows) + 1]] <- data.frame(
    comparison = comp, n_unique_sig_pathways = length(sig_keys))
}

if (length(de_rows) > 0 && length(gse_rows) > 0) {
  de_counts      <- do.call(rbind, de_rows)
  pathway_counts <- do.call(rbind, gse_rows)
  de_counts$comparison      <- factor(de_counts$comparison,      levels = comp_order)
  pathway_counts$comparison <- factor(pathway_counts$comparison, levels = comp_order)
  fwrite(de_counts,      file.path(base_out_dir, "summary_DE_counts.csv"))
  fwrite(pathway_counts, file.path(base_out_dir, "summary_pathway_counts.csv"))
  p_de <- ggplot(de_counts, aes(comparison, n_genes_FDR025, fill = comparison)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n_genes_FDR025), vjust = -0.4, size = 4) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "DE genes (FDR < 0.25)", x = NULL, y = "n genes") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          panel.grid.major.x = element_blank())
  p_gsea <- ggplot(pathway_counts, aes(comparison, n_unique_sig_pathways, fill = comparison)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n_unique_sig_pathways), vjust = -0.4, size = 4) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Unique pathways (padj < 0.05)",
         subtitle = "union across Hallmark + GO_BP × both rankings",
         x = NULL, y = "n pathways") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          panel.grid.major.x = element_blank())
  p_summary <- (p_de | p_gsea) +
    plot_annotation(
      title    = sprintf("Summary — %s", target_label),
      subtitle = sprintf("Clusters: %s", paste(target_clusters, collapse = ",")))
  ggsave(file.path(base_out_dir, "summary.png"), p_summary,
         width = 12, height = 5, dpi = 150)
  cat("\nSummary plots written to:", file.path(base_out_dir, "summary.png"), "\n")
}

cat("\nDone. Outputs in:", base_out_dir, "\n")
