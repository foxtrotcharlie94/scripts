# ============================================================
# GSEA cross-comparison diagnostics + (optional) NES reliability.
#
# Reference comparison: LentiNeg_HBvsLB
# Compared against:      LentiPos_HBvsLB, HB_PosVsNeg, LB_PosVsNeg
#
# PART A (always): from the saved GSEA CSVs in pop_dir, across
#   rankings (logFC, logFCxP) x collection_sets (Hallmark, Hallmark+GO_BP):
#     1. NES scatter (+ Pearson & Spearman)
#     2. Leading-edge Jaccard for pathways significant in both
#     3. Meta-GSEA: ref's significant pathways tested for enrichment in
#        each other comparison's pathway-level NES ranking
#     4. Fraction captured: uncentered NES cos2 (shared fraction),
#        Deming slope (exchange rate), residual (lenti-specific part)
#   -> NES_correlations.csv, metaGSEA_results.csv, captured_fraction.csv,
#      NES_scatter_*.png, captured_summary.png
#
# PART B (optional; run_reliability = TRUE): re-fits edgeR on data splits
#   to estimate each contrast's NES reliability (Spearman-Brown corrected),
#   then rescales the observed cos2 into "fraction of *reliable* signal
#   captured" = cos2 / (rel_ref * rel_other).
#   -> reliability.csv, captured_fraction_ceiling_corrected.csv,
#      captured_ceiling_corrected.png
#   NEEDS THE COUNTS, so set target_clusters to match the population.
#
# Sign conventions (so correlations are interpretable):
#   *_HBvsLB    : positive NES = up in HB
#   *_PosVsNeg  : positive NES = up in LentiPos
# A positive NES correlation between an HBvsLB and a PosVsNeg comparison
# therefore means "HB resembles LentiPos".
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(fgsea)
  library(ggrepel); library(patchwork)
  library(edgeR); library(Matrix); library(msigdbr)
})
n_cores <- 10
data.table::setDTthreads(n_cores)
RhpcBLASctl::blas_set_num_threads(n_cores)
RhpcBLASctl::omp_set_num_threads(n_cores)

# ══════════════════════════════════════════════════════════════════════════════
# Settings
# ══════════════════════════════════════════════════════════════════════════════
pop_dir <- "C:/Users/fc809/Downloads/pseudobulk_DE_GSEA/ClassicalMonocytes_PosVsNegUnpaired"  # population folder

rankings <- c("logFC", "logFCxP")
collection_sets <- list(Hallmark = c("Hallmark"), Hallmark_GO_BP = c("Hallmark","GO_BP"))
ref_comp    <- "LentiNeg_HBvsLB"
other_comps <- c("LentiPos_HBvsLB", "HB_PosVsNeg", "LB_PosVsNeg")
sig_padj    <- 0.05

# GSEA / meta-GSEA params
gsea_species    <- "Mus musculus"
gsea_minSize    <- 15
gsea_maxSize    <- 500
meta_minSize    <- 5
meta_fallback_n <- 15
meta_nPerm      <- 10000
meta_criteria <- list(   # how the reference query set is defined; each runs as its own pass
  padj0.05 = list(col = "padj", cut = 0.05, label = "padj<0.05"),
  padj0.25 = list(col = "padj", cut = 0.25, label = "padj<0.25"),
  pval0.05 = list(col = "pval", cut = 0.05, label = "pval<0.05")
)

# ---- Part B (reliability) ----
run_reliability   <- TRUE
target_clusters   <- c("0")    # MUST match the clusters that generated pop_dir
in_mat_path       <- "C:/Users/fc809/Downloads/count_matrix.txt"
in_meta_path      <- "C:/Users/fc809/Downloads/sample_and_clusterNumber_for_each_cell_res0.5.txt"
n_splits          <- 10        # random splits per contrast (10-20 typical)
split_unit        <- "cell"    # "cell" (technical, conservative) or "donor" (adds biological var)
reliability_nPerm <- 1000      # fewer perms OK: only NES needed
rel_seed          <- 7         # fixed fgsea seed so permutation noise cancels across halves
lenti_gene_name   <- "LentiAll"
THR_POS <- 5e-5; THR_NEG <- 1e-5
min_pct_per_sample_DE <- 0.05
min_cells_per_pb      <- 10
batch_map <- c("HB1"="B1","HB2"="B1","LB1"="B1","LB2"="B1",
               "HB3"="B2","HB4"="B2","LB3"="B2","LB4"="B2","HB5"="B3")

out_dir <- file.path(pop_dir, "_cross_comparison")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════
gsea_path <- function(comp, collection, ranking)
  file.path(pop_dir, comp, sprintf("GSEA_%s_%s_rankBy_%s.csv", collection, comp, ranking))

load_gsea_pooled <- function(comp, collections, ranking) {
  parts <- lapply(collections, function(coll) {
    fp <- gsea_path(comp, coll, ranking)
    if (!file.exists(fp)) { warning("Missing: ", fp); return(NULL) }
    d <- fread(fp); d[, leadingEdge := as.character(leadingEdge)]; d[, collection := coll]; d
  })
  parts <- parts[!sapply(parts, is.null)]
  if (length(parts) == 0) return(NULL)
  rbindlist(parts, fill = TRUE)
}

le_split <- function(x) if (is.na(x) || x == "") character(0) else strsplit(x, ";", fixed = TRUE)[[1]]

make_ref_sets <- function(ref_dt, crit) {
  v   <- ref_dt[[crit$col]]
  sig <- !is.na(v) & v < crit$cut
  up   <- ref_dt[which(sig & ref_dt$NES > 0)][order(-NES), pathway]
  down <- ref_dt[which(sig & ref_dt$NES < 0)][order(NES),  pathway]
  used_fallback <- FALSE
  if (length(up)   < meta_minSize) { up   <- head(ref_dt[NES > 0][order(-NES), pathway], meta_fallback_n); used_fallback <- TRUE }
  if (length(down) < meta_minSize) { down <- head(ref_dt[NES < 0][order(NES),  pathway], meta_fallback_n); used_fallback <- TRUE }
  list(sets = list(ref_up = up, ref_down = down), used_fallback = used_fallback)
}

# Weighted GSEA running enrichment score of a pathway-set within a ranked vector.
# stats_sorted: named numeric (pathway NES of the OTHER comparison), sorted decreasing.
# set_members : pathway names forming the reference query set. Returns running ES per rank.
running_es <- function(stats_sorted, set_members, p = 1) {
  N <- length(stats_sorted)
  inset <- names(stats_sorted) %in% set_members
  ns <- sum(inset)
  if (ns == 0 || ns == N) return(NULL)
  w  <- abs(stats_sorted)^p
  Nr <- sum(w[inset]); if (Nr == 0) return(NULL)
  hit  <- ifelse(inset, w / Nr, 0)
  miss <- ifelse(inset, 0, 1 / (N - ns))
  list(es = cumsum(hit - miss), inset = inset)
}

# Fraction captured between two aligned NES vectors (uncentered: NES=0 is origin)
cap_one <- function(a, b, label) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]; n <- length(a)
  if (n < 5)
    return(data.table(subset = label, n = n, cos2 = NA_real_, signed_cos = NA_real_,
                      slope = NA_real_, residual = NA_real_, pearson = NA_real_))
  Sxx <- sum(a^2); Syy <- sum(b^2); Sxy <- sum(a*b)
  cos_sim <- Sxy / sqrt(Sxx * Syy)
  slope   <- (Syy - Sxx + sqrt((Syy - Sxx)^2 + 4*Sxy^2)) / (2 * Sxy)  # TLS through origin
  data.table(subset = label, n = n, cos2 = cos_sim^2, signed_cos = cos_sim,
             slope = slope, residual = 1 - cos_sim^2, pearson = cor(a, b))
}

# ---- helpers used only by Part B ----
get_msig <- function(species, cat, subcat) {
  a_new <- if (is.null(subcat)) list(species=species, collection=cat) else list(species=species, collection=cat, subcollection=subcat)
  a_old <- if (is.null(subcat)) list(species=species, category=cat)   else list(species=species, category=cat,   subcategory=subcat)
  res <- tryCatch(do.call(msigdbr, a_new), error=function(e) NULL)
  if (is.null(res)) res <- do.call(msigdbr, a_old)
  sc <- intersect(c("gene_symbol","ensembl_gene"), colnames(res))[1]; if (is.na(sc)) sc <- "gene_symbol"
  split(res[[sc]], res$gs_name)
}
build_pb <- function(counts, group_vec, features = NULL) {
  if (!is.null(features)) counts <- counts[features, , drop = FALSE]
  groups <- unique(group_vec)
  M <- sparseMatrix(i = seq_along(group_vec), j = match(group_vec, groups), x = 1,
                    dims = c(length(group_vec), length(groups)), dimnames = list(NULL, groups))
  pb <- as.matrix(counts %*% M); storage.mode(pb) <- "integer"; pb
}
detect_pct_per_sample <- function(counts, sample_vec, min_pct) {
  samples <- unique(sample_vec)
  det <- sapply(samples, function(s) {
    cols <- which(sample_vec == s); if (!length(cols)) return(rep(FALSE, nrow(counts)))
    as.numeric(Matrix::rowSums(counts[, cols, drop = FALSE] > 0)) / length(cols) >= min_pct
  })
  rownames(counts)[rowSums(det) == length(samples)]
}
ucos <- function(x, y) {
  cn <- intersect(names(x), names(y)); x <- x[cn]; y <- y[cn]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5) return(NA_real_)
  sum(x*y) / sqrt(sum(x^2) * sum(y^2))
}
spearman_brown <- function(r) ifelse(is.na(r), NA_real_, ifelse(r <= 0, 0, 2*r/(1+r)))
gsea_nes <- function(ranks, pw) {
  ranks <- ranks[is.finite(ranks)]; ranks <- sort(ranks, decreasing = TRUE)
  if (length(ranks) < gsea_minSize) return(setNames(numeric(0), character(0)))
  set.seed(rel_seed)
  fr <- suppressWarnings(fgsea(pathways = pw, stats = ranks, minSize = gsea_minSize,
                               maxSize = gsea_maxSize, nPermSimple = reliability_nPerm))
  setNames(fr$NES, fr$pathway)
}
rank_from_cells <- function(cnt, md, design) {
  tryCatch({
    gk <- detect_pct_per_sample(cnt, md$sample, min_pct_per_sample_DE)
    if (length(gk) < 50) return(NULL)
    if (design == "paired") {
      md$pb_group <- paste(md$sample, md$genotype, sep = "__")
      gs <- table(md$pb_group); md <- md[md$pb_group %in% names(gs[gs >= min_cells_per_pb]), , drop = FALSE]
      pb <- build_pb(cnt[, rownames(md), drop = FALSE], md$pb_group, features = gk)
      cm <- data.frame(col = colnames(pb), sample = sub("__.*$","",colnames(pb)),
                       genotype = sub("^.*__","",colnames(pb)), stringsAsFactors = FALSE)
      tab <- table(cm$sample, cm$genotype); pd <- rownames(tab)[rowSums(tab > 0) == 2]
      if (length(pd) < 2) return(NULL)
      cm <- cm[cm$sample %in% pd, ]; pb <- pb[, cm$col, drop = FALSE]
      cm$sample <- factor(cm$sample); cm$genotype <- factor(cm$genotype, levels = c("LentiNeg","LentiPos"))
      des <- model.matrix(~ sample + genotype, data = cm); coef <- "genotypeLentiPos"
    } else {
      ss <- table(md$sample); md <- md[md$sample %in% names(ss[ss >= min_cells_per_pb]), , drop = FALSE]
      pb <- build_pb(cnt[, rownames(md), drop = FALSE], md$sample, features = gk)
      cm <- data.frame(col = colnames(pb), sample = colnames(pb), stringsAsFactors = FALSE)
      cm$condition <- ifelse(grepl("^HB", cm$sample), "HB", "LB"); cm$batch <- batch_map[cm$sample]
      if (sum(cm$condition=="HB") < 2 || sum(cm$condition=="LB") < 2) return(NULL)
      cm$condition <- factor(cm$condition, levels = c("LB","HB")); cm$batch <- factor(cm$batch)
      des <- if (nlevels(cm$batch) > 1) model.matrix(~ batch + condition, data = cm) else model.matrix(~ condition, data = cm)
      coef <- "conditionHB"
    }
    y <- DGEList(counts = pb); k <- filterByExpr(y, design = des)
    y <- y[k, , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
    y <- estimateDisp(y, des); fit <- glmQLFit(y, des)
    tt <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
    pP <- pmax(tt$PValue, .Machine$double.xmin)
    list(logFC = setNames(tt$logFC, rownames(tt)),
         logFCxP = setNames(tt$logFC * -log10(pP), rownames(tt)))
  }, error = function(e) NULL)
}
split_halves <- function(md, unit) {
  if (unit == "cell") {
    grp <- paste(md$sample, md$genotype, sep = "__"); a <- logical(nrow(md))
    for (g in unique(grp)) { idx <- which(grp == g); s <- sample(idx); a[s[seq_len(floor(length(s)/2))]] <- TRUE }
  } else {
    sa <- character()
    for (cond in unique(md$condition)) { cs <- sample(unique(md$sample[md$condition == cond])); sa <- c(sa, cs[seq_len(floor(length(cs)/2))]) }
    a <- md$sample %in% sa
  }
  list(A = a, B = !a)
}

# ══════════════════════════════════════════════════════════════════════════════
# PART A — cross-comparison from saved GSEA CSVs
# ══════════════════════════════════════════════════════════════════════════════
all_corr <- list(); all_meta <- list(); all_captured <- list(); all_meta_curves <- list()

for (ranking in rankings) {
  for (set_name in names(collection_sets)) {
    colls <- collection_sets[[set_name]]
    cat(sprintf("\n================ %s | rankBy %s ================\n", set_name, ranking))

    ref_dt <- load_gsea_pooled(ref_comp, colls, ranking)
    if (is.null(ref_dt)) { cat("  (reference missing — skipped)\n"); next }
    ref_sets_by_crit <- lapply(meta_criteria, function(cr) make_ref_sets(ref_dt, cr))
    for (cn2 in names(meta_criteria)) {
      s <- ref_sets_by_crit[[cn2]]
      cat(sprintf("  ref sets @%s: ref_up=%d  ref_down=%d%s\n", meta_criteria[[cn2]]$label,
                  length(s$sets$ref_up), length(s$sets$ref_down),
                  if (s$used_fallback) "  (fallback to top/bottom N)" else ""))
    }

    scatter_plots <- list()
    for (oc in other_comps) {
      other_dt <- load_gsea_pooled(oc, colls, ranking)
      if (is.null(other_dt)) next

      # 1. NES scatter + correlations
      m <- merge(ref_dt[,   .(pathway, collection, NES_ref = NES, padj_ref = padj)],
                 other_dt[, .(pathway, NES_oth = NES, padj_oth = padj)], by = "pathway")
      m <- m[is.finite(NES_ref) & is.finite(NES_oth)]
      pear  <- if (nrow(m) > 2) cor(m$NES_ref, m$NES_oth) else NA_real_
      spear <- if (nrow(m) > 2) cor(m$NES_ref, m$NES_oth, method = "spearman") else NA_real_
      all_corr[[length(all_corr)+1]] <- data.frame(
        ranking = ranking, collection_set = set_name, ref = ref_comp, other = oc,
        n_pathways = nrow(m), pearson = pear, spearman = spear)

      # 4. Fraction captured (all pathways + union-significant)
      m_sig <- m[padj_ref < sig_padj | padj_oth < sig_padj]
      cap <- rbind(cap_one(m$NES_ref, m$NES_oth, "all_pathways"),
                   cap_one(m_sig$NES_ref, m_sig$NES_oth, "union_significant"))
      cap[, `:=`(ranking = ranking, collection_set = set_name, ref = ref_comp, other = oc)]
      all_captured[[length(all_captured)+1]] <- cap

      m[, sig_cat := fifelse(padj_ref < sig_padj & padj_oth < sig_padj, "both",
                      fifelse(padj_ref < sig_padj, "ref only",
                       fifelse(padj_oth < sig_padj, "other only", "ns")))]
      lab <- head(m[padj_ref < sig_padj & padj_oth < sig_padj][order(-(abs(NES_ref)+abs(NES_oth)))], 12)
      rng <- range(c(m$NES_ref, m$NES_oth), na.rm = TRUE)
      scatter_plots[[oc]] <- ggplot(m, aes(NES_ref, NES_oth)) +
        geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(aes(color = sig_cat, shape = collection), size = 1.8, alpha = 0.7) +
        geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
        geom_text_repel(data = lab, aes(label = pathway), size = 2.4, max.overlaps = 20, segment.color = "grey70") +
        scale_color_manual(values = c("both"="firebrick","ref only"="steelblue","other only"="darkorange","ns"="grey75"), name = NULL) +
        scale_shape_manual(values = c("Hallmark"=17,"GO_BP"=16), name = NULL) +
        coord_equal(xlim = rng, ylim = rng) +
        labs(title = sprintf("NES  %s  vs  %s", ref_comp, oc),
             subtitle = sprintf("Pearson r = %.2f | Spearman = %.2f | n = %d", pear, spear, nrow(m)),
             x = sprintf("NES — %s", ref_comp), y = sprintf("NES — %s", oc)) +
        theme_bw(base_size = 10)

      # 2. Leading-edge Jaccard
      common_sig <- m[padj_ref < sig_padj & padj_oth < sig_padj, pathway]
      if (length(common_sig) > 0) {
        je <- rbindlist(lapply(common_sig, function(pw) {
          le_r <- le_split(ref_dt[pathway == pw, leadingEdge][1]); le_o <- le_split(other_dt[pathway == pw, leadingEdge][1])
          inter <- length(intersect(le_r, le_o)); uni <- length(union(le_r, le_o))
          data.table(pathway = pw, NES_ref = ref_dt[pathway == pw, NES][1], NES_oth = other_dt[pathway == pw, NES][1],
                     n_le_ref = length(le_r), n_le_oth = length(le_o), n_shared = inter,
                     jaccard = if (uni>0) inter/uni else NA_real_,
                     shared_genes = paste(sort(intersect(le_r, le_o)), collapse = ";"))
        }))[order(-jaccard)]
        fwrite(je, file.path(out_dir, sprintf("leadingEdge_jaccard_%s_rankBy_%s_%s_vs_%s.csv", set_name, ranking, ref_comp, oc)))
        cat(sprintf("  %s vs %s: %d pathways sig in both | median LE Jaccard = %.2f\n",
                    ref_comp, oc, nrow(je), median(je$jaccard, na.rm = TRUE)))
      } else cat(sprintf("  %s vs %s: no pathways sig in both\n", ref_comp, oc))

      # 3. Meta-GSEA — once per inclusion criterion; accumulate running-ES curves
      stats <- other_dt$NES; names(stats) <- other_dt$pathway
      stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
      for (crit_name in names(meta_criteria)) {
        crit_lab <- meta_criteria[[crit_name]]$label
        ref_sets <- ref_sets_by_crit[[crit_name]]$sets
        sets_use <- ref_sets[sapply(ref_sets, function(s) sum(s %in% names(stats)) >= meta_minSize)]
        if (length(sets_use) == 0 || length(stats) < meta_minSize) {
          cat(sprintf("  meta-GSEA %s vs %s @%s: insufficient set/stat overlap\n", ref_comp, oc, crit_lab)); next
        }
        set.seed(42)
        mg <- as.data.table(fgsea(pathways = sets_use, stats = stats, minSize = meta_minSize,
                                  maxSize = length(stats), nPermSimple = meta_nPerm))
        mg[, `:=`(ranking = ranking, collection_set = set_name, criterion = crit_lab, ref = ref_comp, other = oc)]
        mg[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
        all_meta[[length(all_meta)+1]] <- mg
        cat(sprintf("  meta-GSEA %s vs %s @%s:\n", ref_comp, oc, crit_lab))
        for (i in seq_len(nrow(mg)))
          cat(sprintf("    %-9s NES=%+.2f  padj=%.2e  (set size %d)\n", mg$pathway[i], mg$NES[i], mg$padj[i], mg$size[i]))

        # Running-ES curves for plotting (one per direction set present)
        for (dir_name in names(sets_use)) {
          r <- running_es(stats, sets_use[[dir_name]])
          if (is.null(r)) next
          nes_d  <- mg[pathway == dir_name, NES]; padj_d <- mg[pathway == dir_name, padj]
          all_meta_curves[[length(all_meta_curves)+1]] <- data.table(
            ranking = ranking, collection_set = set_name, criterion = crit_lab, other = oc,
            direction = dir_name, rank = seq_along(stats), es = r$es, inset = r$inset,
            nes = if (length(nes_d)) nes_d else NA_real_,
            padj = if (length(padj_d)) padj_d else NA_real_)
        }
      }
    }

    if (length(scatter_plots) > 0) {
      fig <- wrap_plots(scatter_plots, ncol = length(scatter_plots)) +
        plot_annotation(title = sprintf("NES cross-comparison — %s", basename(pop_dir)),
                        subtitle = sprintf("collection set: %s | rankBy: %s", set_name, ranking))
      ggsave(file.path(out_dir, sprintf("NES_scatter_%s_rankBy_%s.png", set_name, ranking)),
             fig, width = 6 * length(scatter_plots), height = 6.3, dpi = 150)
    }
  }
}

if (length(all_corr) > 0) {
  corr_dt <- rbindlist(lapply(all_corr, as.data.table))
  fwrite(corr_dt, file.path(out_dir, "NES_correlations.csv"))
  cat("\nNES correlation summary:\n"); print(corr_dt)
}
if (length(all_meta) > 0) {
  meta_dt <- rbindlist(all_meta, fill = TRUE)
  setcolorder(meta_dt, c("ranking","collection_set","criterion","ref","other","pathway","NES","pval","padj","size"))
  fwrite(meta_dt, file.path(out_dir, "metaGSEA_results.csv"))
}

# ── Meta-GSEA running-score enrichment plots ──────────────────────────────────
if (length(all_meta_curves) > 0) {
  curves <- rbindlist(all_meta_curves, fill = TRUE)
  curves[, other := factor(other, levels = other_comps)]
  dir_cols <- c(ref_up = "firebrick", ref_down = "steelblue")

  # annotation label per panel (NES / padj for each direction)
  anno <- unique(curves[, .(ranking, collection_set, criterion, other, direction, nes, padj)])
  anno[, lab := sprintf("%s: NES=%+.2f, padj=%.2g", sub("ref_","",direction), nes, padj)]

  # (a) Detailed file per collection_set x ranking x criterion: facet direction ~ other, with rug
  for (sn in names(collection_sets)) for (rk in rankings) for (cn3 in names(meta_criteria)) {
    crit_lab <- meta_criteria[[cn3]]$label
    d <- curves[collection_set == sn & ranking == rk & criterion == crit_lab]
    if (nrow(d) == 0) next
    a <- anno[collection_set == sn & ranking == rk & criterion == crit_lab]
    p <- ggplot(d, aes(rank, es)) +
      geom_hline(yintercept = 0, color = "grey70") +
      geom_line(aes(color = direction), linewidth = 0.6) +
      geom_rug(data = d[inset == TRUE], aes(color = direction), sides = "b", alpha = 0.5) +
      geom_text(data = a, aes(label = lab, color = direction),
                x = Inf, y = Inf, hjust = 1.02, vjust = 1.4, size = 4.8, fontface = "bold", show.legend = FALSE) +
      facet_grid(direction ~ other, scales = "free_y") +
      scale_color_manual(values = dir_cols, name = NULL) +
      labs(title = sprintf("Meta-GSEA enrichment — %s | %s | rankBy %s | ref set %s",
                           basename(pop_dir), sn, rk, crit_lab),
           subtitle = sprintf("ref = %s | x = each comparison's pathway rank, ticks = reference set members", ref_comp),
           x = "rank of pathways in the other comparison (NES, high \u2192 low)", y = "running ES") +
      theme_bw(base_size = 10) +
      theme(legend.text = element_text(size = 15),
            legend.key.size = unit(1.4, "lines"),
            strip.text = element_text(size = 12))
    ggsave(file.path(out_dir, sprintf("metaGSEA_enrichment_%s_rankBy_%s_%s.png", sn, rk, cn3)),
           p, width = 13, height = 6, dpi = 150)
  }

  # (b) One summary contact sheet: rows = set|rank|criterion, cols = other, both directions overlaid
  curves[, config := factor(sprintf("%s | %s | %s", collection_set, ranking, criterion))]
  anno[,   config := factor(sprintf("%s | %s | %s", collection_set, ranking, criterion),
                            levels = levels(curves$config))]
  anno[, ypos := Inf]
  anno[, vj   := ifelse(direction == "ref_up", 1.3, 2.9)]
  ps <- ggplot(curves, aes(rank, es, color = direction)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.45) +
    geom_text(data = anno, aes(label = sprintf("%s %+.2f", sub("ref_","",direction), nes), y = ypos, vjust = vj),
              x = Inf, hjust = 1.03, size = 2.1, show.legend = FALSE) +
    facet_grid(config ~ other, scales = "free_x") +
    scale_color_manual(values = dir_cols, name = NULL) +
    labs(title = sprintf("Meta-GSEA enrichment summary — %s", basename(pop_dir)),
         subtitle = sprintf("ref = %s | each panel: reference up/down pathway-set enrichment within the other comparison's pathway ranking",
                            ref_comp),
         x = "rank of pathways in the other comparison (NES, high \u2192 low)", y = "running ES") +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0, size = 7), legend.position = "top")
  ggsave(file.path(out_dir, "metaGSEA_enrichment_summary.png"), ps,
         width = 12, height = 2 + 1.5 * length(levels(curves$config)), dpi = 150, limitsize = FALSE)
  cat("\nMeta-GSEA enrichment plots written.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# PART A2 — gene-level meta-GSEA
# Reference DE signature (significant up/down GENES) tested for enrichment within
# each other comparison's ranked GENE list. Same idea as the pathway meta-GSEA,
# but members are genes (from DE_*.csv) so there is no collection dimension.
# ══════════════════════════════════════════════════════════════════════════════
gene_criteria <- list(   # how the reference gene signature is defined
  FDR0.05  = list(col = "FDR",    cut = 0.05, label = "FDR<0.05"),
  FDR0.25  = list(col = "FDR",    cut = 0.25, label = "FDR<0.25"),
  pval0.05 = list(col = "PValue", cut = 0.05, label = "PValue<0.05")
)
gene_rank_col     <- c(logFC = "logFC", logFCxP = "logFC_x_negLogP")  # other comparison's gene ranking
gene_meta_minSize <- 10

load_de <- function(comp) {
  fp <- file.path(pop_dir, comp, sprintf("DE_%s_byPValue.csv", comp))
  if (!file.exists(fp)) { warning("Missing DE: ", fp); return(NULL) }
  fread(fp)
}
make_gene_sets <- function(de, crit) {
  v <- de[[crit$col]]; sig <- !is.na(v) & v < crit$cut
  list(ref_up   = de[which(sig & de$logFC > 0)][order(-logFC), gene],
       ref_down = de[which(sig & de$logFC < 0)][order(logFC),  gene])
}

all_gene_meta <- list(); all_gene_curves <- list()
ref_de <- load_de(ref_comp)
if (!is.null(ref_de)) {
  cat("\n================ Gene-level meta-GSEA ================\n")
  for (rk in rankings) {
    rcol <- gene_rank_col[[rk]]
    gene_sets_by_crit <- lapply(gene_criteria, function(cr) make_gene_sets(ref_de, cr))
    for (cn in names(gene_criteria))
      cat(sprintf("  [rankBy %s] ref genes @%s: up=%d down=%d\n", rk, gene_criteria[[cn]]$label,
                  length(gene_sets_by_crit[[cn]]$ref_up), length(gene_sets_by_crit[[cn]]$ref_down)))
    for (oc in other_comps) {
      other_de <- load_de(oc); if (is.null(other_de)) next
      stats <- other_de[[rcol]]; names(stats) <- other_de$gene
      stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
      for (cn in names(gene_criteria)) {
        crit_lab <- gene_criteria[[cn]]$label
        gs <- gene_sets_by_crit[[cn]]
        sets_use <- gs[sapply(gs, function(s) sum(s %in% names(stats)) >= gene_meta_minSize)]
        if (length(sets_use) == 0 || length(stats) < gene_meta_minSize) {
          cat(sprintf("  gene meta %s vs %s [rankBy %s] @%s: insufficient overlap\n", ref_comp, oc, rk, crit_lab)); next
        }
        set.seed(42)
        mg <- as.data.table(fgsea(pathways = sets_use, stats = stats,
                                  minSize = gene_meta_minSize, maxSize = length(stats),
                                  nPermSimple = meta_nPerm))
        mg[, `:=`(ranking = rk, criterion = crit_lab, ref = ref_comp, other = oc)]
        mg[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
        all_gene_meta[[length(all_gene_meta)+1]] <- mg
        cat(sprintf("  gene meta %s vs %s [rankBy %s] @%s:\n", ref_comp, oc, rk, crit_lab))
        for (i in seq_len(nrow(mg)))
          cat(sprintf("    %-9s NES=%+.2f  padj=%.2e  (set size %d)\n", mg$pathway[i], mg$NES[i], mg$padj[i], mg$size[i]))
        for (dir_name in names(sets_use)) {
          r <- running_es(stats, sets_use[[dir_name]]); if (is.null(r)) next
          nes_d <- mg[pathway == dir_name, NES]; padj_d <- mg[pathway == dir_name, padj]
          all_gene_curves[[length(all_gene_curves)+1]] <- data.table(
            ranking = rk, criterion = crit_lab, other = oc, direction = dir_name,
            rank = seq_along(stats), es = r$es, inset = r$inset,
            nes = if (length(nes_d)) nes_d else NA_real_,
            padj = if (length(padj_d)) padj_d else NA_real_)
        }
      }
    }
  }
}

if (length(all_gene_meta) > 0) {
  gmeta_dt <- rbindlist(all_gene_meta, fill = TRUE)
  setcolorder(gmeta_dt, c("ranking","criterion","ref","other","pathway","NES","pval","padj","size"))
  fwrite(gmeta_dt, file.path(out_dir, "geneMetaGSEA_results.csv"))

  gcurves <- rbindlist(all_gene_curves, fill = TRUE)
  gcurves[, other := factor(other, levels = other_comps)]
  dir_cols <- c(ref_up = "firebrick", ref_down = "steelblue")
  ganno <- unique(gcurves[, .(ranking, criterion, other, direction, nes, padj)])
  ganno[, lab := sprintf("%s: NES=%+.2f, padj=%.2g", sub("ref_","",direction), nes, padj)]

  # thin long gene curves for rendering (keep peak + endpoints)
  thin_curve <- function(d, maxpts = 2500) {
    if (nrow(d) <= maxpts) return(d)
    keep <- sort(unique(c(1L, nrow(d), which.max(abs(d$es)),
                          as.integer(round(seq(1, nrow(d), length.out = maxpts))))))
    d[keep]
  }
  gcurves_thin <- gcurves[, thin_curve(.SD), by = .(ranking, criterion, other, direction)]

  # (a) detailed file per ranking x criterion: facet direction ~ other
  for (rk in rankings) for (cn in names(gene_criteria)) {
    crit_lab <- gene_criteria[[cn]]$label
    d <- gcurves_thin[ranking == rk & criterion == crit_lab]; if (nrow(d) == 0) next
    a <- ganno[ranking == rk & criterion == crit_lab]
    drug <- gcurves[ranking == rk & criterion == crit_lab & inset == TRUE]
    p <- ggplot(d, aes(rank, es)) +
      geom_hline(yintercept = 0, color = "grey70") +
      geom_line(aes(color = direction), linewidth = 0.6) +
      geom_rug(data = drug, aes(color = direction), sides = "b", alpha = 0.25, linewidth = 0.2) +
      geom_text(data = a, aes(label = lab, color = direction),
                x = Inf, y = Inf, hjust = 1.02, vjust = 1.4, size = 2.8, show.legend = FALSE) +
      facet_grid(direction ~ other, scales = "free_y") +
      scale_color_manual(values = dir_cols, name = NULL) +
      labs(title = sprintf("Gene-level meta-GSEA — %s | rankBy %s | ref signature %s",
                           basename(pop_dir), rk, crit_lab),
           subtitle = sprintf("ref = %s | x = each comparison's gene rank (high \u2192 low)", ref_comp),
           x = "rank of genes in the other comparison", y = "running ES") +
      theme_bw(base_size = 10)
    ggsave(file.path(out_dir, sprintf("geneMetaGSEA_enrichment_rankBy_%s_%s.png", rk, cn)),
           p, width = 13, height = 6, dpi = 150)
  }

  # (b) one summary contact sheet: rows = ranking | criterion, cols = other
  gcurves_thin[, config := factor(sprintf("%s | %s", ranking, criterion))]
  ganno[, config := factor(sprintf("%s | %s", ranking, criterion), levels = levels(gcurves_thin$config))]
  ganno[, ypos := Inf]
  ganno[, vj   := ifelse(direction == "ref_up", 1.3, 2.9)]
  gps <- ggplot(gcurves_thin, aes(rank, es, color = direction)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.45) +
    geom_text(data = ganno, aes(label = sprintf("%s %+.2f", sub("ref_","",direction), nes), y = ypos, vjust = vj),
              x = Inf, hjust = 1.03, size = 2.1, show.legend = FALSE) +
    facet_grid(config ~ other, scales = "free_x") +
    scale_color_manual(values = dir_cols, name = NULL) +
    labs(title = sprintf("Gene-level meta-GSEA enrichment summary — %s", basename(pop_dir)),
         subtitle = sprintf("ref = %s | each panel: reference up/down GENE signature within the other comparison's gene ranking", ref_comp),
         x = "rank of genes in the other comparison (high \u2192 low)", y = "running ES") +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0, size = 7), legend.position = "top")
  ggsave(file.path(out_dir, "geneMetaGSEA_enrichment_summary.png"), gps,
         width = 12, height = 2 + 1.5 * length(levels(gcurves_thin$config)), dpi = 150, limitsize = FALSE)
  cat("Gene-level meta-GSEA plots written.\n")
}

cap_dt <- data.table()
if (length(all_captured) > 0) {
  cap_dt <- rbindlist(all_captured, fill = TRUE)
  setcolorder(cap_dt, c("ranking","collection_set","ref","other","subset","n","cos2","residual","signed_cos","slope","pearson"))
  fwrite(cap_dt, file.path(out_dir, "captured_fraction.csv"))

  axis_type <- c(LentiPos_HBvsLB = "within-axis (burden)", HB_PosVsNeg = "cross-axis (lenti)", LB_PosVsNeg = "cross-axis (lenti)")
  axis_cols <- c("cross-axis (lenti)" = "firebrick", "within-axis (burden)" = "grey70")
  pdat <- cap_dt[subset == "all_pathways"]
  pdat[, other := factor(other, levels = other_comps)]
  pdat[, axis := factor(axis_type[as.character(other)], levels = names(axis_cols))]

  pA <- ggplot(pdat, aes(other, cos2, fill = axis)) +
    geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.2f", cos2)), vjust = -0.3, size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    facet_grid(collection_set ~ ranking) +
    scale_fill_manual(values = axis_cols, name = NULL) +
    scale_y_continuous(limits = c(0, 1.08), expand = expansion(mult = c(0, 0))) +
    labs(title = "Shared pathway-signature fraction (NES cos\u00b2)",
         subtitle = sprintf("ref = %s | headroom to 1.0 = residual (not on burden axis)", ref_comp),
         x = NULL, y = "captured fraction (cos\u00b2)") +
    theme_bw(base_size = 11) + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))
  pB <- ggplot(pdat, aes(other, slope, color = axis)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey60") +
    geom_segment(aes(xend = other, yend = 0), linewidth = 0.5) + geom_point(size = 3) +
    geom_text(aes(label = sprintf("%.2f", slope)), vjust = -0.7, size = 3) +
    facet_grid(collection_set ~ ranking) +
    scale_color_manual(values = axis_cols, name = NULL) +
    labs(title = "Deming slope (exchange rate): NES_other per unit NES_ref",
         subtitle = "sign = direction (>0 concordant) | dotted = 1:1 | dashed = 0", x = NULL, y = "Deming slope") +
    theme_bw(base_size = 11) + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(out_dir, "captured_summary.png"), pA / pB +
         plot_annotation(title = sprintf("How much of the lenti signature lies on the burden axis — %s", basename(pop_dir))),
         width = 11, height = 11, dpi = 150)
}

# ══════════════════════════════════════════════════════════════════════════════
# PART B — split-half NES reliability + noise-ceiling correction (optional)
# ══════════════════════════════════════════════════════════════════════════════
if (run_reliability && nrow(cap_dt) > 0) {
  cat("\n================ Reliability (Part B) ================\n")

  cat("Loading cluster assignments...\n")
  meta_cl <- fread(in_meta_path)
  setnames(meta_cl, c("cells","harmony_clusters_res0.5"), c("cell","cluster"))
  meta_cl[, `:=`(cell = as.character(cell), cluster = as.character(cluster))]
  cells_target <- meta_cl[cluster %in% target_clusters, cell]

  cat("Reading count_matrix (selected cells)...\n")
  hdr <- names(fread(in_mat_path, nrows = 0)); gene_col <- hdr[1]
  cells_to_read <- intersect(cells_target, hdr[-1])
  dt <- fread(in_mat_path, select = c(gene_col, cells_to_read))
  genes <- as.character(dt[[1]]); mat <- as.matrix(dt[, -1, with = FALSE]); rownames(mat) <- genes
  counts <- as(mat, "CsparseMatrix"); rm(dt, mat); gc()

  lenti_umis <- as.numeric(counts[lenti_gene_name, ]); total_umis <- as.numeric(Matrix::colSums(counts))
  ratio <- lenti_umis / pmax(total_umis, 1)
  genotype <- ifelse(ratio >= THR_POS, "LentiPos", ifelse(ratio <= THR_NEG, "LentiNeg", NA))
  counts <- counts[rownames(counts) != lenti_gene_name, , drop = FALSE]
  keepc <- !is.na(genotype); counts <- counts[, keepc, drop = FALSE]
  meta <- data.frame(row.names = colnames(counts), cell = colnames(counts), stringsAsFactors = FALSE)
  meta$sample <- sub("_.*$", "", meta$cell)
  meta$condition <- ifelse(grepl("^HB", meta$sample), "HB", "LB")
  meta$batch <- batch_map[meta$sample]; meta$genotype <- genotype[keepc]
  cat(sprintf("Loaded %d genes x %d classified cells\n", nrow(counts), ncol(counts)))

  collection_pw <- list(Hallmark = get_msig(gsea_species, "H", NULL),
                        GO_BP    = get_msig(gsea_species, "C5", "GO:BP"))
  contrasts_spec <- list(
    LentiNeg_HBvsLB = list(cells = function(md) md$genotype == "LentiNeg", design = "unpaired"),
    LentiPos_HBvsLB = list(cells = function(md) md$genotype == "LentiPos", design = "unpaired"),
    HB_PosVsNeg     = list(cells = function(md) md$condition == "HB",      design = "paired"),
    LB_PosVsNeg     = list(cells = function(md) md$condition == "LB",      design = "paired"))

  rel_rows <- list()
  for (cn in names(contrasts_spec)) {
    spec <- contrasts_spec[[cn]]; idx <- spec$cells(meta)
    cnt_c <- counts[, idx, drop = FALSE]; md_c <- meta[idx, , drop = FALSE]
    cat(sprintf("\n[%s] %d cells | %d %s-level splits\n", cn, ncol(cnt_c), n_splits, split_unit))
    acc <- list()
    for (sp in seq_len(n_splits)) {
      set.seed(1000 + sp); h <- split_halves(md_c, split_unit)
      rA <- rank_from_cells(cnt_c[, h$A, drop=FALSE], md_c[h$A,,drop=FALSE], spec$design)
      rB <- rank_from_cells(cnt_c[, h$B, drop=FALSE], md_c[h$B,,drop=FALSE], spec$design)
      if (is.null(rA) || is.null(rB)) { cat(sprintf("  split %d skipped\n", sp)); next }
      for (mt in rankings) {
        nesA <- lapply(collection_pw, function(pw) gsea_nes(rA[[mt]], pw))
        nesB <- lapply(collection_pw, function(pw) gsea_nes(rB[[mt]], pw))
        for (sn in names(collection_sets)) {
          cols <- collection_sets[[sn]]
          acc[[paste(mt, sn, sep="|")]] <- c(acc[[paste(mt, sn, sep="|")]],
                                             ucos(do.call(c, nesA[cols]), do.call(c, nesB[cols])))
        }
      }
    }
    for (key in names(acc)) {
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]; rh <- mean(acc[[key]], na.rm = TRUE)
      rel_rows[[length(rel_rows)+1]] <- data.table(
        contrast = cn, ranking = parts[1], collection_set = parts[2],
        n_splits_used = sum(!is.na(acc[[key]])), r_half_mean = rh,
        r_half_sd = sd(acc[[key]], na.rm = TRUE), reliability = spearman_brown(rh))
    }
  }
  rel_dt <- rbindlist(rel_rows)
  fwrite(rel_dt, file.path(out_dir, "reliability.csv"))
  cat("\nReliability (Spearman-Brown):\n")
  print(rel_dt[, .(contrast, ranking, collection_set, reliability = round(reliability,3))])

  cap <- copy(cap_dt[subset == "all_pathways"])
  cap <- merge(cap, rel_dt[, .(ref = contrast, ranking, collection_set, rel_ref = reliability)],
               by = c("ref","ranking","collection_set"), all.x = TRUE, sort = FALSE)
  cap <- merge(cap, rel_dt[, .(other = contrast, ranking, collection_set, rel_other = reliability)],
               by = c("other","ranking","collection_set"), all.x = TRUE, sort = FALSE)
  cap[, ceiling_cos2 := rel_ref * rel_other]
  cap[, captured_of_reliable := pmin(cos2 / ceiling_cos2, 1)]
  cap[, corrected_signed_cor := pmax(pmin(signed_cos / sqrt(ceiling_cos2), 1), -1)]
  fwrite(cap, file.path(out_dir, "captured_fraction_ceiling_corrected.csv"))
  cat("\nCeiling-corrected captured fraction:\n")
  print(cap[, .(ranking, collection_set, other, observed_cos2 = round(cos2,3),
                ceiling = round(ceiling_cos2,3), captured_of_reliable = round(captured_of_reliable,3))])

  axis_type <- c(LentiPos_HBvsLB = "within-axis (burden)", HB_PosVsNeg = "cross-axis (lenti)", LB_PosVsNeg = "cross-axis (lenti)")
  axis_cols <- c("cross-axis (lenti)" = "firebrick", "within-axis (burden)" = "grey70")
  cap[, other := factor(other, levels = other_comps)]
  cap[, axis := factor(axis_type[as.character(other)], levels = names(axis_cols))]
  p <- ggplot(cap, aes(other)) +
    geom_col(aes(y = cos2, fill = axis), width = 0.7) +
    geom_errorbar(aes(ymin = ceiling_cos2, ymax = ceiling_cos2), width = 0.7, color = "black", linewidth = 0.6) +
    geom_text(aes(y = cos2, label = sprintf("%.0f%%", 100*captured_of_reliable)), vjust = -0.4, size = 3) +
    facet_grid(collection_set ~ ranking) +
    scale_fill_manual(values = axis_cols, name = NULL) +
    scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0,0))) +
    labs(title = sprintf("Observed shared fraction vs noise ceiling — %s", basename(pop_dir)),
         subtitle = paste0("bar = observed cos\u00b2 | black line = ceiling (rel_ref \u00d7 rel_other) | label = % of *reliable* signal captured\n",
                           sprintf("ref = %s | %d %s-level splits", ref_comp, n_splits, split_unit)),
         x = NULL, y = "cos\u00b2") +
    theme_bw(base_size = 11) + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(out_dir, "captured_ceiling_corrected.png"), p, width = 11, height = 8, dpi = 150)
}

cat("\nDone. Outputs in:", out_dir, "\n")
