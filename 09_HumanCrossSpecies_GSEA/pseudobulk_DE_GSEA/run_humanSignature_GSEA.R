#!/usr/bin/env Rscript
# =============================================================================
# Human-signature GSEA in mouse monocyte DGE
# -----------------------------------------------------------------------------
# Gene sets  : from each of 4 human DGE files, among genes with padj/fdr<cutoff,
#              the TOP_N most up- and most down-regulated by logFC  -> UP/DOWN sets.
# Ranked list: mouse pseudobulk edgeR DGE, one per (monocyte subset x comparison),
#              mapped mouse->human orthologs (babelgene), ranked two ways:
#                metric 1  "logFC"          : logFC
#                metric 2  "logFC_signedP"  : logFC * -log10(FDR)
#
# Asks: are the human TET2/CHIP monocyte UP/DOWN programs coordinately
#       enriched in each mouse monocyte contrast?  (positive NES on an UP set
#       = that human-up program is also up in the mouse ranking.)
# =============================================================================

# ============================== CONFIG =======================================
HUMAN_FILES <- c(
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/DESeq2_Cohort_CCUS_CHIP_vs_Control/CD14+_Monocytes_TET2_CHIP_vs_Control.tsv",
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/DESeq2_Cohort_CH_vs_Control/CD14+_Monocytes_TET2_vs_none.tsv",
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/MAST_DE_Mutant_vs_Control/TET2 CD14 Mutant vs Control.csv",
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/MAST_DE_Mutant_vs_Wildtype/TET2 CD14 Mutant vs Wildtype.csv"
)

MOUSE_ROOT  <- "C:/Users/fc809/Downloads/pseudobulk_DE_GSEA"
SUBSETS     <- c("ClassicalMonocytes", "NonClassicalMonocytes", "MHCII+Monocytes")
COMPARISONS <- c("LentiNeg_HBvsLB", "LentiPos_HBvsLB", "HB_PosVsNeg", "LB_PosVsNeg")

OUTPUT_DIR  <- "C:/Users/fc809/Downloads/pseudobulk_DE_GSEA/_humanSignature_GSEA"

SIG_CUTOFF  <- 0.05     # human padj/fdr cutoff for signature membership
TOP_N       <- 200      # per signature: take the TOP_N most up- and most down-regulated
                        #   (by logFC) among genes passing SIG_CUTOFF
MIN_SIZE    <- 5
MAX_SIZE    <- 5000
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("fgsea", quietly = TRUE))        BiocManager::install("fgsea", update = FALSE, ask = FALSE)
  if (!requireNamespace("babelgene", quietly = TRUE))    install.packages("babelgene", repos = "https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE))   install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE))      install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("patchwork", quietly = TRUE))    install.packages("patchwork", repos = "https://cloud.r-project.org")
  library(fgsea); library(babelgene); library(data.table); library(ggplot2); library(patchwork)
})
FIG_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------- flexible column detection ---------------------------
pick <- function(cn, candidates) {
  hit <- candidates[tolower(candidates) %in% tolower(cn)]
  if (length(hit)) cn[match(tolower(hit[1]), tolower(cn))] else NA_character_
}
read_dge <- function(path) {
  ext <- tolower(tools::file_ext(path))
  d <- fread(path, sep = if (ext == "tsv") "\t" else ",", data.table = FALSE)
  cn <- names(d)
  gcol <- pick(cn, c("gene","Gene","symbol","SYMBOL","primerid","gene_symbol","names","X","V1"))
  if (is.na(gcol) && (cn[1] %in% c("V1","X",""))) gcol <- cn[1]
  lcol <- pick(cn, c("log2FoldChange","avg_log2FC","logFC","coef","avg_logFC","log2FC"))
  acol <- pick(cn, c("padj","p_val_adj","FDR","fdr","adj.P.Val","qvalue"))
  pcol <- pick(cn, c("PValue","pvalue","p_val","P.Value","Pr(>Chisq)","pval"))
  if (is.na(gcol) || is.na(lcol))
    stop(sprintf("No gene/logFC column in %s | cols: %s", basename(path), paste(cn, collapse=", ")))
  padj <- if (!is.na(acol)) d[[acol]] else if (!is.na(pcol)) d[[pcol]] else NA
  data.frame(gene = as.character(d[[gcol]]), logFC = as.numeric(d[[lcol]]),
             padj = as.numeric(padj), stringsAsFactors = FALSE)
}

# ------------------ build UP / DOWN human signatures -------------------------
message("Building human UP/DOWN signatures (padj<", SIG_CUTOFF, ", top/bottom ", TOP_N, " by logFC) ...")
signatures <- list()
for (path in HUMAN_FILES) {
  d   <- read_dge(path)
  tag <- gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(path)))
  sig <- d[!is.na(d$padj) & d$padj < SIG_CUTOFF & !is.na(d$logFC), ]
  sig <- sig[!duplicated(sig$gene), ]
  up  <- sig[sig$logFC > 0, ]; up <- up[order(-up$logFC), ]          # most up-regulated first
  dn  <- sig[sig$logFC < 0, ]; dn <- dn[order( dn$logFC), ]          # most down-regulated first
  signatures[[paste0(tag, "__UP")]]   <- head(up$gene, TOP_N)
  signatures[[paste0(tag, "__DOWN")]] <- head(dn$gene, TOP_N)
  message(sprintf("  %-50s UP=%4d DOWN=%4d", substr(tag,1,50),
                  length(signatures[[paste0(tag,"__UP")]]),
                  length(signatures[[paste0(tag,"__DOWN")]])))
}
signatures <- signatures[lengths(signatures) >= MIN_SIZE]

# ---------------- mouse->human ortholog table (build once) -------------------
# (we map per ranked list, but caching the full table speeds things up)
mouse_human <- NULL
get_map <- function(genes) {
  o <- babelgene::orthologs(genes = genes, species = "mouse", human = FALSE)
  o[, c("symbol", "human_symbol")]            # symbol = mouse, human_symbol = human
}

mouse_to_human_ranks <- function(d, metric) {
  d <- d[!is.na(d$gene) & !is.na(d$logFC), , drop = FALSE]
  if (metric == "logFC") {
    d$val <- d$logFC
  } else {
    p <- d$padj; p[is.na(p)] <- 1
    nz <- suppressWarnings(min(p[p > 0])); if (is.finite(nz)) p[p == 0] <- nz
    d$val <- d$logFC * -log10(p)
  }
  d <- d[is.finite(d$val), ]
  map <- get_map(unique(d$gene))
  m <- merge(d[, c("gene","val")], map, by.x = "gene", by.y = "symbol")
  m <- m[!is.na(m$human_symbol) & nzchar(m$human_symbol) & is.finite(m$val), ]   # drop unmapped
  if (!nrow(m)) return(setNames(numeric(0), character(0)))
  m <- m[order(-abs(m$val)), ]
  m <- m[!duplicated(m$human_symbol), ]      # collapse paralogs -> strongest signal
  v <- sort(setNames(m$val, m$human_symbol), decreasing = TRUE)
  v[!is.na(names(v)) & nzchar(names(v))]     # final guard: no NA/empty names for fgsea
}

# locate the per-comparison DE file (any of the 3 sorted variants = same data)
find_de <- function(subset, comp) {
  dir <- file.path(MOUSE_ROOT, subset, comp)
  hits <- list.files(dir, pattern = sprintf("^DE_%s_.*\\.csv$", comp), full.names = TRUE)
  if (!length(hits)) return(NA_character_)
  pref <- hits[grepl("byPValue", hits)]
  if (length(pref)) pref[1] else hits[1]
}

# --------------------------------- run ---------------------------------------
metrics <- c("logFC", "logFC_signedP")
all_res     <- list()
ranks_store <- list()   # keyed "subset||comp||metric" -> ranked vector (for plots)

for (subset in SUBSETS) for (comp in COMPARISONS) {
  de <- find_de(subset, comp)
  if (is.na(de)) { warning(sprintf("Missing DE file: %s / %s", subset, comp)); next }
  d <- read_dge(de)
  for (metric in metrics) {
    ranks <- mouse_to_human_ranks(d, metric)
    message(sprintf("%-22s | %-16s | %-13s | %d genes", subset, comp, metric, length(ranks)))
    if (length(ranks) < MIN_SIZE) { warning("  too few mapped genes; skipping"); next }
    res <- fgsea(pathways = signatures, stats = ranks,
                 minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0)
    res <- res[order(res$padj), ]
    res$leadingEdge  <- vapply(res$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    res$subset <- subset; res$comparison <- comp; res$metric <- metric
    fwrite(res, file.path(OUTPUT_DIR,
           sprintf("%s__%s__%s.csv", gsub("[^A-Za-z0-9]+","",subset), comp, metric)))
    all_res[[length(all_res)+1]] <- res
    ranks_store[[paste(subset, comp, metric, sep = "||")]] <- ranks
  }
}

if (!length(all_res)) stop("No results - check the configured paths.")

combined <- rbindlist(all_res, use.names = TRUE, fill = TRUE)
setcolorder(combined, c("subset","comparison","metric","pathway","pval","padj","NES","ES","size"))
fwrite(combined, file.path(OUTPUT_DIR, "humanSignature_GSEA_all_results.csv"))
message(sprintf("\nGSEA done. %d signatures x %d runs. Now drawing figures ...",
                length(signatures), length(all_res)))

# =============================================================================
#                               VISUALS
# =============================================================================
# colors matching the reference: DOWN = red, UP = blue
COL_DOWN <- "#C0392B"; COL_UP <- "#1F3A93"

# human file short labels (drop the __UP/__DOWN suffix)
file_tags <- unique(sub("__(UP|DOWN)$", "", names(signatures)))
short_lab <- function(tag) {
  s <- gsub("_+", " ", tag); if (nchar(s) > 38) paste0(substr(s,1,36),"..") else s
}

# ---- classic weighted running enrichment score (gseaParam = 1) --------------
running_es <- function(stats, pathway) {
  stats <- sort(stats, decreasing = TRUE)
  hits  <- names(stats) %in% pathway
  Nh <- sum(hits); n <- length(stats)
  if (Nh == 0) return(NULL)
  r <- abs(stats)
  Phit  <- cumsum(ifelse(hits, r, 0)) / sum(r[hits])
  Pmiss <- cumsum(ifelse(!hits, 1, 0)) / (n - Nh)
  list(curve = data.frame(rank = seq_len(n), ES = Phit - Pmiss),
       ticks = which(hits), metric = as.numeric(stats))
}

# ---- 3-panel enrichment plot: UP + DOWN of one file on one ranked list -------
enrichment_plot <- function(ranks, upSet, downSet, lab, nesU, fdrU, nesD, fdrD, ttl) {
  eU <- running_es(ranks, upSet); eD <- running_es(ranks, downSet)
  if (is.null(eU) || is.null(eD)) return(NULL)
  n <- length(ranks)

  curve <- rbind(data.frame(eU$curve, grp = "UP"), data.frame(eD$curve, grp = "DOWN"))
  legU <- sprintf("%s-up  NES=%.2f FDR=%.2g", lab, nesU, fdrU)
  legD <- sprintf("%s-down  NES=%.2f FDR=%.2g", lab, nesD, fdrD)
  curve$grp <- factor(curve$grp, levels = c("DOWN","UP"), labels = c(legD, legU))

  p_es <- ggplot(curve, aes(rank, ES, color = grp)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = .3) +
    geom_line(linewidth = .7) +
    scale_color_manual(values = setNames(c(COL_DOWN, COL_UP), c(legD, legU))) +
    labs(y = "Running Enrichment Score", color = NULL, title = ttl) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 11, face = "plain"))

  ticks <- rbind(data.frame(rank = eD$ticks, y0 = 0.55, y1 = 1.0, col = COL_DOWN),
                 data.frame(rank = eU$ticks, y0 = 0.0,  y1 = 0.45, col = COL_UP))
  p_tk <- ggplot(ticks) +
    geom_segment(aes(x = rank, xend = rank, y = y0, yend = y1), color = ticks$col, linewidth = .25) +
    scale_x_continuous(limits = c(1, n)) + scale_y_continuous(limits = c(0,1)) +
    theme_void()

  wf <- data.frame(rank = seq_len(n), metric = eD$metric)
  p_wf <- ggplot(wf, aes(rank, metric)) +
    geom_area(fill = "grey75") +
    labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
    theme_classic(base_size = 11)

  p_es / p_tk / p_wf + plot_layout(heights = c(3, 0.9, 1.4))
}

# ---- one multipage PDF per human file (page = subset x comparison x metric) --
for (tag in file_tags) {
  upN <- paste0(tag, "__UP"); dnN <- paste0(tag, "__DOWN")
  if (is.null(signatures[[upN]]) || is.null(signatures[[dnN]])) next
  pdf(file.path(FIG_DIR, sprintf("enrichment__%s.pdf", gsub("[^A-Za-z0-9]+","_",tag))),
      width = 8.5, height = 7.5, onefile = TRUE)
  for (key in names(ranks_store)) {
    kk <- strsplit(key, "\\|\\|")[[1]]; .subset <- kk[1]; .comp <- kk[2]; .metric <- kk[3]
    ranks <- ranks_store[[key]]
    # compute the row filter OUTSIDE data.table's [ ] so subset/metric refer to the
    # loop values, not the identically-named columns (data.table i-scoping gotcha)
    mask <- combined$pathway %in% c(upN, dnN) & combined$subset == .subset &
            combined$comparison == .comp & combined$metric == .metric
    r <- combined[mask, ]
    gv <- function(p, col) { v <- r[[col]][r$pathway==p]; if (length(v)) v[1] else NA }
    pl <- enrichment_plot(ranks, signatures[[upN]], signatures[[dnN]], short_lab(tag),
                          gv(upN,"NES"), gv(upN,"padj"), gv(dnN,"NES"), gv(dnN,"padj"),
                          sprintf("%s: %s (%s) [%s]", short_lab(tag), .comp, .subset, .metric))
    if (!is.null(pl)) print(pl)
  }
  dev.off()
}
message("  enrichment PDFs written.")

# ---- helper labels for summary plots ----------------------------------------
combined$context <- paste(combined$subset, combined$comparison, sep = " | ")
combined$sig     <- combined$pathway
combined$star    <- ifelse(combined$padj < 0.05, "*", "")

# ==== (1) NES heatmap : signatures x contexts, faceted by metric =============
hm <- ggplot(combined, aes(context, sig, fill = NES)) +
  geom_tile(color = "white", linewidth = .3) +
  geom_text(aes(label = star), vjust = .78, size = 5) +
  facet_wrap(~ metric, ncol = 1) +
  scale_fill_gradient2(low = COL_UP, mid = "white", high = COL_DOWN, midpoint = 0) +
  labs(x = NULL, y = NULL, title = "Human-signature NES across mouse monocyte contrasts",
       subtitle = "* padj < 0.05") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "summary_1_NES_heatmap.pdf"), hm, width = 11, height = 9)

# ==== (2) UP vs DOWN directional concordance =================================
combined$file <- sub("__(UP|DOWN)$", "", combined$pathway)
combined$dir  <- ifelse(grepl("__UP$", combined$pathway), "UP", "DOWN")
wide <- dcast(combined, subset + comparison + metric + file ~ dir, value.var = "NES")
ud <- ggplot(wide, aes(UP, DOWN, color = metric)) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_abline(slope = -1, intercept = 0, linetype = 2, color = "grey60") +
  geom_point(size = 2.4, alpha = .8) +
  labs(title = "Directional concordance: UP vs DOWN signature NES",
       subtitle = "Coherent signal falls near the dashed anti-diagonal (UP up <-> DOWN down)",
       x = "NES (UP signature)", y = "NES (DOWN signature)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG_DIR, "summary_2_UPvsDOWN.pdf"), ud, width = 7.5, height = 6.5)

# ==== (3) cross-context NES correlation (per metric) =========================
for (mt in unique(combined$metric)) {
  m <- as.data.frame(dcast(combined[combined$metric == mt, ], sig ~ context, value.var = "NES"))
  M <- as.matrix(m[, -1, drop = FALSE]); rownames(M) <- m$sig
  cc <- cor(M, use = "pairwise.complete.obs")
  cl <- as.data.frame(as.table(cc)); names(cl) <- c("c1","c2","r")
  g <- ggplot(cl, aes(c1, c2, fill = r)) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.6) +
    scale_fill_gradient2(low = COL_UP, mid = "white", high = COL_DOWN, midpoint = 0, limits = c(-1,1)) +
    labs(title = sprintf("Cross-context NES correlation [%s]", mt), x = NULL, y = NULL) +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(FIG_DIR, sprintf("summary_3_context_corr_%s.pdf", mt)), g, width = 8.5, height = 7.5)
}

# ==== (4) metric agreement: NES(logFC) vs NES(logFC_signedP) =================
ma <- dcast(combined, subset + comparison + sig ~ metric, value.var = "NES")
if (all(c("logFC","logFC_signedP") %in% names(ma))) {
  rval <- suppressWarnings(cor(ma$logFC, ma$logFC_signedP, use = "complete.obs"))
  g4 <- ggplot(ma, aes(logFC, logFC_signedP)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
    geom_point(aes(color = subset), size = 2.2, alpha = .8) +
    labs(title = "Metric agreement of NES",
         subtitle = sprintf("Pearson r = %.3f", rval),
         x = "NES (rank by logFC)", y = "NES (rank by logFC x -log10FDR)") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG_DIR, "summary_4_metric_agreement.pdf"), g4, width = 7.5, height = 6.5)
}

message(sprintf("Done. Results + figures in %s/  (figures/ holds enrichment PDFs + 4 summaries).", OUTPUT_DIR))
