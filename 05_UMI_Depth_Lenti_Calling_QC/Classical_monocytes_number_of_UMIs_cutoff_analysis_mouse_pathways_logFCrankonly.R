## ============================================================
##  GSEA sensitivity analysis — Classical monocytes
##  4 comparisons x 4 UMI-filtering cases x 2 FDR thresholds
##  Output: 4x2 grid of Hallmark heatmaps
## ============================================================

library(data.table)
library(tidyverse)
library(edgeR)
library(fgsea)
library(msigdbr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

setwd("C:/Users/fc809/Downloads/")

## ── Parameters ────────────────────────────────────────────────────────────────
MATRIX_FILE  <- "classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt"
THRESH_FILE  <- "Lenti_detection_thresholds_SCTv2.csv"
LENTI_FREQ   <- 4e-5      # LentiPos: Lenti UMIs / total UMIs >= this
FIXED_THRESH  <- 13700    # fixed nUMI threshold (cases 3 & 4)
FIXED_THRESH2 <- 9227     # pooled NB-derived threshold (cases 5 & 6)
MIN_CELLS    <- 100       # min cells per sample per group
FDR_LO       <- 0.25     # lenient FDR (left panels)
FDR_HI       <- 0.1      # strict FDR  (right panels)
N_PERM       <- 5000
FDR_CUT      <- 0.05     # FDR for master ComplexHeatmap

## ── Load Hallmark gene sets ───────────────────────────────────────────────────
message("Loading Hallmark gene sets...")
gs_hallmark <- msigdbr(species = "Mus musculus", collection = "MH", db_species = "MM") %>%
  { split(.$gene_symbol, .$gs_name) } %>%
  lapply(unique)
message("  ", length(gs_hallmark), " gene sets loaded")

## ── Load thresholds ───────────────────────────────────────────────────────────
thresholds <- read.csv(THRESH_FILE)
thresh_vec <- setNames(thresholds$UMI_Threshold, thresholds$Sample)

## ── Parse gene x cell matrix ──────────────────────────────────────────────────
message("\nReading matrix...")
mat_raw    <- as.data.frame(fread(MATRIX_FILE))
gene_names <- mat_raw[[1]]
mat        <- as.matrix(mat_raw[, -1])
rownames(mat) <- gene_names
barcodes   <- colnames(mat)
message("  ", nrow(mat), " genes x ", ncol(mat), " cells")

## ── Per-cell metadata ─────────────────────────────────────────────────────────
lenti_row <- which(grepl("^[Ll]enti$", gene_names))[1]
if (is.na(lenti_row)) lenti_row <- which(grepl("[Ll]enti", gene_names))[1]
message("  Lenti row: '", gene_names[lenti_row], "'")

numi_vec    <- colSums(mat)
lenti_umis  <- mat[lenti_row, ]
lenti_freq  <- lenti_umis / numi_vec
sample_vec  <- sub("_.*", "", barcodes)
cond_vec    <- ifelse(grepl("^HB", sample_vec), "HB", "LB")

meta_all <- data.frame(
  barcode      = barcodes,
  sample       = sample_vec,
  condition    = cond_vec,
  nUMI         = numi_vec,
  lenti_freq   = lenti_freq,
  Lenti_status = ifelse(lenti_freq >= LENTI_FREQ, "LentiPos", "LentiNeg"),
  stringsAsFactors = FALSE
)
message("  HB cells: ", sum(meta_all$condition == "HB"),
        " | LB cells: ", sum(meta_all$condition == "LB"))
message("  LentiPos: ", sum(meta_all$Lenti_status == "LentiPos"),
        " | LentiNeg: ", sum(meta_all$Lenti_status == "LentiNeg"))

## ── Define 4 comparisons ──────────────────────────────────────────────────────
# Each comparison specifies:
#   cell_filter : which rows of meta_all to keep before case filtering
#   group_var   : column used as the contrast variable
#   ref_level   : reference level (denominator)
#   test_level  : test level (numerator)
#   paired      : include sample in edgeR design (TRUE for within-sample comparisons)
#   label       : plot label

comparisons <- list(
  list(
    id         = "LentiPosNeg_LB",
    label      = "LentiPos vs LentiNeg (LB)",
    cell_filter = function(m) m[m$condition == "LB", ],
    group_var   = "Lenti_status",
    ref_level   = "LentiNeg",
    test_level  = "LentiPos",
    paired      = TRUE,   # same samples have both groups → block by sample
    do_match    = TRUE    # UMI matching makes sense (LentiPos has higher nUMI)
  ),
  list(
    id         = "LentiPosNeg_HB",
    label      = "LentiPos vs LentiNeg (HB)",
    cell_filter = function(m) m[m$condition == "HB", ],
    group_var   = "Lenti_status",
    ref_level   = "LentiNeg",
    test_level  = "LentiPos",
    paired      = TRUE,
    do_match    = TRUE
  ),
  list(
    id         = "HBvsLB_LentiNeg",
    label      = "HB vs LB (LentiNeg)",
    cell_filter = function(m) m[m$Lenti_status == "LentiNeg", ],
    group_var   = "condition",
    ref_level   = "LB",
    test_level  = "HB",
    paired      = FALSE,
    do_match    = TRUE
  ),
  list(
    id         = "HBvsLB_LentiPos",
    label      = "HB vs LB (LentiPos)",
    cell_filter = function(m) m[m$Lenti_status == "LentiPos", ],
    group_var   = "condition",
    ref_level   = "LB",
    test_level  = "HB",
    paired      = FALSE,
    do_match    = TRUE
  )
)

## ── Define 4 UMI-filtering cases ─────────────────────────────────────────────
cases <- list(
  "1_all_cells" = list(
    label      = "Case 1: all cells",
    umi_filter = "none",
    umi_match  = FALSE
  ),
  "2_all_matched" = list(
    label      = "Case 2: all cells + UMI matched",
    umi_filter = "none",
    umi_match  = TRUE
  ),
  "3_fixed_thresh" = list(
    label      = paste0("Case 3: nUMI \u2265 ", FIXED_THRESH),
    umi_filter = "fixed",
    umi_match  = FALSE
  ),
  "4_fixed_matched" = list(
    label      = paste0("Case 4: nUMI \u2265 ", FIXED_THRESH, " + UMI matched"),
    umi_filter = "fixed",
    umi_match  = TRUE
  ),
  "5_pooled_thresh" = list(
    label      = paste0("Case 5: nUMI \u2265 ", FIXED_THRESH2),
    umi_filter = "fixed2",
    umi_match  = FALSE
  ),
  "6_pooled_matched" = list(
    label      = paste0("Case 6: nUMI \u2265 ", FIXED_THRESH2, " + UMI matched"),
    umi_filter = "fixed2",
    umi_match  = TRUE
  ),
  "7_sample_thresh" = list(
    label      = "Case 7: sample-specific SCTv2",
    umi_filter = "sample",
    umi_match  = FALSE
  ),
  "8_sample_matched" = list(
    label      = "Case 8: sample-specific SCTv2 + UMI matched",
    umi_filter = "sample",
    umi_match  = TRUE
  )
)

## ── Helper: filter samples with < MIN_CELLS ───────────────────────────────────
filter_min_cells <- function(meta_sub, group_var, paired) {
  counts <- meta_sub %>%
    group_by(sample, .data[[group_var]]) %>%
    summarise(n = n(), .groups = "drop")
  
  if (paired) {
    wide <- counts %>%
      pivot_wider(names_from = all_of(group_var), values_from = n, values_fill = 0)
    lvls <- unique(meta_sub[[group_var]])
    keep <- wide
    for (lv in lvls) keep <- keep[keep[[lv]] >= MIN_CELLS, ]
    keep_samples <- keep$sample
  } else {
    keep_samples <- counts %>%
      dplyr::filter(n >= MIN_CELLS) %>%
      pull(sample) %>% unique()
  }
  
  message("    Samples surviving filter: ", paste(sort(keep_samples), collapse = ", "))
  meta_sub[meta_sub$sample %in% keep_samples, ]
}

## ── Helper: UMI matching (subsample test group DOWN to ref distribution) ──────
match_umi <- function(meta_sub, group_var, ref_level, test_level) {
  result <- list()
  for (samp in unique(meta_sub$sample)) {
    s    <- meta_sub[meta_sub$sample == samp, ]
    ref  <- s[s[[group_var]] == ref_level,  ]
    test <- s[s[[group_var]] == test_level, ]
    if (nrow(ref) == 0 || nrow(test) == 0) { result[[samp]] <- s; next }
    breaks <- unique(quantile(ref$nUMI, probs = seq(0, 1, 0.1)))
    if (length(breaks) < 2) {
      keep_test <- test[sample(nrow(test), min(nrow(test), nrow(ref))), ]
    } else {
      ref$bin  <- cut(ref$nUMI,  breaks = breaks, include.lowest = TRUE)
      test$bin <- cut(test$nUMI, breaks = breaks, include.lowest = TRUE)
      bin_counts <- table(ref$bin)
      keep_test <- do.call(rbind, lapply(names(bin_counts), function(b) {
        avail  <- test[!is.na(test$bin) & test$bin == b, ]
        n_want <- bin_counts[b]
        if (nrow(avail) == 0 || n_want == 0) return(NULL)
        avail[sample(nrow(avail), min(nrow(avail), n_want)), ]
      }))
      ref$bin       <- NULL
      keep_test$bin <- NULL
    }
    result[[samp]] <- rbind(ref, keep_test)
  }
  do.call(rbind, result)
}

## ── Helper: pseudobulk edgeR ──────────────────────────────────────────────────
run_edger <- function(meta_sub, count_mat, group_var, ref_level, test_level, paired) {
  levels_grp <- c(ref_level, test_level)
  pb_list <- list(); pb_meta <- list()
  for (samp in unique(meta_sub$sample)) {
    for (grp in levels_grp) {
      cells <- meta_sub$barcode[meta_sub$sample == samp & meta_sub[[group_var]] == grp]
      if (length(cells) == 0) next
      cn <- paste0(samp, "_", grp)
      pb_list[[cn]] <- rowSums(count_mat[, cells, drop = FALSE])
      pb_meta[[cn]] <- data.frame(sample = samp, group = grp, stringsAsFactors = FALSE)
    }
  }
  if (length(pb_list) == 0) return(NULL)
  pb_counts <- do.call(cbind, pb_list)
  pb_df     <- do.call(rbind, pb_meta)
  rownames(pb_df) <- colnames(pb_counts)
  pb_df$group  <- factor(pb_df$group,  levels = levels_grp)
  pb_df$sample <- factor(pb_df$sample)
  n_ref  <- sum(pb_df$group == ref_level)
  n_test <- sum(pb_df$group == test_level)
  message("    Pseudobulk: ", n_ref, " ", ref_level, " + ", n_test, " ", test_level)
  if (n_ref < 2 || n_test < 2) {
    message("    WARNING: <2 replicates — skipping"); return(NULL)
  }
  design <- if (paired) model.matrix(~ sample + group, data = pb_df) else
    model.matrix(~ group,          data = pb_df)
  coef_name <- paste0("group", test_level)
  if (!coef_name %in% colnames(design)) {
    message("    WARNING: coefficient '", coef_name, "' not found"); return(NULL)
  }
  dge <- DGEList(counts = pb_counts)
  dge <- dge[filterByExpr(dge, design), ]
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  res <- glmQLFTest(fit, coef = coef_name)
  tt  <- topTags(res, n = Inf)$table
  tt$gene <- rownames(tt)
  tt
}

## ── Helper: fgsea ─────────────────────────────────────────────────────────────
run_fgsea <- function(deg) {
  if (is.null(deg)) return(NULL)
  ranked <- setNames(deg$logFC, deg$gene)
  ranked <- ranked[!duplicated(names(ranked)) & is.finite(ranked)]
  set.seed(42)
  set.seed(42)
  fgsea(pathways = gs_hallmark, stats = ranked,
        minSize = 15, maxSize = 500, nPermSimple = N_PERM, nproc = 1) %>%
    as_tibble() %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
}

## ── Helper: build one heatmap ─────────────────────────────────────────────────
make_hm <- function(nm, fdr_label, title_str) {
  if (is.null(nm) || nrow(nm) == 0 || ncol(nm) == 0) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label=paste0("No gene sets\n(FDR < ", fdr_label, ")"),
                      size=12, color="grey50") + theme_void() +
             labs(title=title_str) +
             theme(plot.title=element_text(face="bold",size=18,hjust=0.5)))
  }
  mx <- max(abs(nm), na.rm=TRUE); if (!is.finite(mx)||mx==0) mx <- 1
  hm_df <- as.data.frame(nm) %>%
    rownames_to_column("case") %>%
    pivot_longer(-case, names_to="pathway", values_to="NES") %>%
    mutate(case    = factor(case,    levels=rev(rownames(nm))),
           pathway = factor(pathway, levels=colnames(nm)))
  ggplot(hm_df, aes(x=pathway, y=case, fill=NES)) +
    geom_tile(color="white", linewidth=0.4) +
    scale_fill_gradientn(
      colours  = c("#2980B9","white","#C0392B"),
      values   = scales::rescale(c(-mx, 0, mx)),
      limits   = c(-mx, mx),
      na.value = "grey90",
      name     = "NES"
    ) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    labs(title=title_str,
         subtitle=paste0("grey = FDR \u2265 ", fdr_label),
         x=NULL, y=NULL) +
    theme_minimal(base_size=18) +
    theme(
      plot.title      = element_text(face="bold", size=20, hjust=0.5),
      plot.subtitle   = element_text(size=15, hjust=0.5, color="grey40"),
      axis.text.x     = element_text(angle=45, hjust=1, size=14),
      axis.text.y     = element_text(size=14),
      legend.title    = element_text(size=13),
      legend.text     = element_text(size=12),
      legend.position = "right",
      plot.margin     = margin(t=10, r=15, b=5, l=15)
    )
}

## ── Helper: build NES matrix for a given results df and FDR cutoff ────────────
build_nes_matrix <- function(results_df, case_labels, fdr_cut) {
  sig <- results_df %>% dplyr::filter(padj < fdr_cut) %>%
    pull(pathway) %>% unique()
  if (length(sig) == 0) return(NULL)
  nw <- results_df %>% dplyr::filter(pathway %in% sig) %>%
    dplyr::select(case, pathway, NES) %>%
    pivot_wider(names_from=pathway, values_from=NES) %>%
    column_to_rownames("case")
  pw <- results_df %>% dplyr::filter(pathway %in% sig) %>%
    dplyr::select(case, pathway, padj) %>%
    pivot_wider(names_from=pathway, values_from=padj) %>%
    column_to_rownames("case")
  nw  <- nw[ intersect(case_labels, rownames(nw)),  , drop=FALSE]
  pw  <- pw[ intersect(case_labels, rownames(pw)),  , drop=FALSE]
  for (r in rownames(nw))
    for (c in colnames(nw))
      if (is.na(pw[r,c]) || pw[r,c] >= fdr_cut) nw[r,c] <- NA
  colnames(nw) <- str_to_title(gsub("_"," ",sub("^HALLMARK_","",colnames(nw))))
  tmp <- nw; tmp[is.na(tmp)] <- 0
  if (ncol(tmp)>1) nw <- nw[, hclust(dist(t(tmp)))$order, drop=FALSE]
  as.matrix(nw)
}

## ── Cell count bar plots ─────────────────────────────────────────────────────
message("Building cell count plots...")

cell_counts_all <- list()

for (cmp in comparisons) {
  for (case_id in names(cases)) {
    cfg  <- cases[[case_id]]
    meta <- cmp$cell_filter(meta_all)
    
    if (cfg$umi_filter == "fixed") {
      meta <- meta[meta$nUMI >= FIXED_THRESH, ]
    } else if (cfg$umi_filter == "fixed2") {
      meta <- meta[meta$nUMI >= FIXED_THRESH2, ]
    }
    
    meta <- filter_min_cells(meta, cmp$group_var, cmp$paired)
    if (nrow(meta) == 0) next
    
    if (cfg$umi_match && cmp$do_match) {
      meta <- match_umi(meta, cmp$group_var, cmp$ref_level, cmp$test_level)
      meta <- filter_min_cells(meta, cmp$group_var, cmp$paired)
      if (nrow(meta) == 0) next
    }
    
    counts <- meta %>%
      group_by(sample, .data[[cmp$group_var]]) %>%
      summarise(n_cells = n(), .groups = "drop") %>%
      mutate(
        comparison = cmp$label,
        case       = cfg$label,
        group      = .data[[cmp$group_var]]
      )
    cell_counts_all[[paste(cmp$id, case_id)]] <- counts
  }
}

cell_counts_df <- bind_rows(cell_counts_all)

# Ordered factor levels — must be set BEFORE creating x_label
case_levels  <- sapply(cases, `[[`, "label")
group_levels <- c("LentiPos","LentiNeg","HB","LB")

cell_counts_df <- cell_counts_df %>%
  mutate(
    case       = factor(case,  levels = case_levels),
    group      = factor(group, levels = group_levels),
    comparison = factor(comparison, levels = sapply(comparisons, `[[`, "label"))
  ) %>%
  arrange(comparison, group, case)

# Build x_label with correct order per comparison
# Use numeric position to guarantee ordering within each facet
cell_counts_df <- cell_counts_df %>%
  mutate(
    case_num  = as.integer(case),
    group_num = as.integer(group),
    # short case label: strip "Case N: "
    case_short = gsub("^Case \\d+: ", "", as.character(case)),
    x_label = paste0(as.character(group), "\n", case_short)
  )

# Build global x_label factor with correct ordering
x_order <- cell_counts_df %>%
  distinct(group, group_num, case, case_num, x_label) %>%
  arrange(group_num, case_num) %>%
  pull(x_label) %>%
  unique()

cell_counts_df$x_label <- factor(cell_counts_df$x_label, levels = x_order)

# Sample colour palette
all_samples <- sort(unique(cell_counts_df$sample))
sample_pal  <- setNames(scales::hue_pal()(length(all_samples)), all_samples)

# Number of unique x positions
n_x <- length(x_order)

# Use case (short) on x-axis, facet by comparison × group
cell_counts_df <- cell_counts_df %>%
  mutate(
    case_short2 = paste0("C", as.integer(case)),  # C1..C8 for compact x labels
    facet_label = paste0(comparison, "\n", group)
  )

# Order facet labels: comparison first, then group within comparison
facet_order <- cell_counts_df %>%
  distinct(comparison, group, group_num) %>%
  arrange(comparison, group_num) %>%
  mutate(fl = paste0(comparison, "\n", group)) %>%
  pull(fl)

cell_counts_df$facet_label <- factor(cell_counts_df$facet_label,
                                     levels = unique(facet_order))

n_facet_cols <- length(unique(cell_counts_df$group))  # 2 groups side by side
n_facet_rows <- length(comparisons)

p_cells <- ggplot(cell_counts_df,
                  aes(x = case_short2, y = n_cells, fill = sample)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.15) +
  facet_wrap(~ facet_label, nrow = n_facet_rows, ncol = n_facet_cols,
             scales = "free_y") +
  scale_fill_manual(values = sample_pal, name = "Sample") +
  scale_x_discrete(labels = function(x) {
    # map C1..C8 back to short names
    lbl <- c(C1="all",C2="all\n+match",C3="≥13700",C4="≥13700\n+match",
             C5="≥9227",C6="≥9227\n+match",C7="SCTv2",C8="SCTv2\n+match")
    lbl[x]
  }) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Cells per sample per case — Classical Monocytes",
    subtitle = "Each bar stacked by sample  |  rows = comparisons, columns = groups  |  C1–C8 = cases",
    x = NULL, y = "N cells"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 15),
    plot.subtitle   = element_text(size = 11, color = "grey40"),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 11),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(size = 12),
    strip.text      = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.spacing   = unit(1.0, "lines"),
    plot.margin     = margin(t = 10, r = 10, b = 10, l = 10)
  )

plot_h <- 4 * n_facet_rows + 2
plot_w <- 5 * n_facet_cols + 3

ggsave("GSEA_Monocytes_cell_counts_per_case.pdf",
       p_cells, width = plot_w, height = plot_h, limitsize = FALSE)
ggsave("GSEA_Monocytes_cell_counts_per_case.png",
       p_cells, width = plot_w, height = plot_h, dpi = 150, limitsize = FALSE)

message("Cell count plot saved.")

## ── Main loop: comparisons x cases ───────────────────────────────────────────
case_labels <- sapply(cases, `[[`, "label")
all_plots_lo <- list()   # FDR < 0.25 panels
all_plots_hi <- list()   # FDR < 0.10 panels

for (cmp in comparisons) {
  message("\n\n========================================")
  message("Comparison: ", cmp$label)
  
  gsea_results <- list()
  
  for (case_id in names(cases)) {
    cfg <- cases[[case_id]]
    message("\n  --- ", cfg$label, " ---")
    
    # 1. Filter to relevant cells for this comparison
    meta <- cmp$cell_filter(meta_all)
    
    # 2. nUMI filter
    if (cfg$umi_filter == "fixed") {
      meta <- meta[meta$nUMI >= FIXED_THRESH, ]
      message("  After nUMI >= ", FIXED_THRESH, ": ", nrow(meta), " cells")
    } else if (cfg$umi_filter == "fixed2") {
      meta <- meta[meta$nUMI >= FIXED_THRESH2, ]
      message("  After nUMI >= ", FIXED_THRESH2, ": ", nrow(meta), " cells")
    } else if (cfg$umi_filter == "sample") {
      meta <- meta[meta$nUMI >= thresh_vec[meta$sample], ]
      message("  After sample-specific thresholds: ", nrow(meta), " cells")
    }
    
    # 3. Min cells filter
    meta <- filter_min_cells(meta, cmp$group_var, cmp$paired)
    if (nrow(meta) == 0) { message("  No cells — skipping"); next }
    
    # 4. UMI matching (only if comparison supports it)
    if (cfg$umi_match && cmp$do_match) {
      meta <- match_umi(meta, cmp$group_var, cmp$ref_level, cmp$test_level)
      message("  After UMI matching: ", nrow(meta), " cells")
      meta <- filter_min_cells(meta, cmp$group_var, cmp$paired)
      if (nrow(meta) == 0) { message("  No cells after matching — skipping"); next }
    } else if (cfg$umi_match && !cmp$do_match) {
      message("  (UMI matching skipped for this comparison type)")
    }
    
    # 5. edgeR + fgsea
    deg <- run_edger(meta, mat, cmp$group_var, cmp$ref_level,
                     cmp$test_level, cmp$paired)
    res <- run_fgsea(deg)
    if (!is.null(res)) {
      gsea_results[[case_id]] <- res %>% mutate(case = cfg$label)
    }
  }
  
  if (length(gsea_results) == 0) {
    message("  No results for ", cmp$label)
    all_plots_lo[[cmp$id]] <- make_hm(NULL, FDR_LO, cmp$label)
    all_plots_hi[[cmp$id]] <- make_hm(NULL, FDR_HI, cmp$label)
    next
  }
  
  all_res <- bind_rows(gsea_results)
  write.csv(all_res,
            paste0("GSEA_Monocytes_", cmp$id, "_all_results.csv"),
            row.names = FALSE)
  
  nm_lo <- build_nes_matrix(all_res, case_labels, FDR_LO)
  nm_hi <- build_nes_matrix(all_res, case_labels, FDR_HI)
  
  all_plots_lo[[cmp$id]] <- make_hm(nm_lo, FDR_LO, cmp$label)
  all_plots_hi[[cmp$id]] <- make_hm(nm_hi, FDR_HI, cmp$label)
  
  n_sig_lo <- if (!is.null(nm_lo)) ncol(nm_lo) else 0
  n_sig_hi <- if (!is.null(nm_hi)) ncol(nm_hi) else 0
  message("  Sig gene sets: FDR<", FDR_LO, ": ", n_sig_lo,
          " | FDR<", FDR_HI, ": ", n_sig_hi)
}

## ── Assemble 4x2 grid ─────────────────────────────────────────────────────────
message("\n=== Building combined plot ===")

left_col  <- wrap_plots(lapply(comparisons, function(c) all_plots_lo[[c$id]]),
                        ncol=1) +
  plot_annotation(title = paste0("FDR < ", FDR_LO),
                  theme = theme(plot.title = element_text(face="bold", size=22,
                                                          hjust=0.5)))

right_col <- wrap_plots(lapply(comparisons, function(c) all_plots_hi[[c$id]]),
                        ncol=1) +
  plot_annotation(title = paste0("FDR < ", FDR_HI),
                  theme = theme(plot.title = element_text(face="bold", size=22,
                                                          hjust=0.5)))

combined <- left_col | right_col

ggsave("GSEA_Monocytes_sensitivity_heatmap.pdf",
       combined, width=52, height=28, limitsize=FALSE)
ggsave("GSEA_Monocytes_sensitivity_heatmap.png",
       combined, width=52, height=28, dpi=100, limitsize=FALSE)
## ── Combined heatmap with dendrogram ──────────────────────────────────────────
message("\n=== Building combined heatmap with dendrogram ===")

# Collect all GSEA results across all comparisons
all_results_combined <- list()
for (cmp in comparisons) {
  csv_file <- paste0("GSEA_Monocytes_", cmp$id, "_all_results.csv")
  if (file.exists(csv_file)) {
    res <- read.csv(csv_file, stringsAsFactors = FALSE)
    res$comparison <- cmp$label
    all_results_combined[[cmp$id]] <- res
  }
}

if (length(all_results_combined) > 0) {
  all_combined <- bind_rows(all_results_combined)
  
  # Build combined NES matrix for each FDR threshold
  build_combined_heatmap <- function(fdr_cut, fdr_label) {
    
    # Row ID = "comparison | case label" — keeps blocks contiguous
    all_combined <- all_combined %>%
      mutate(row_id = paste0(comparison, "  |  ", case))
    
    # Significant pathways across ALL comparisons at this FDR
    sig_paths <- all_combined %>%
      dplyr::filter(padj < fdr_cut) %>%
      pull(pathway) %>% unique()
    
    if (length(sig_paths) == 0) {
      message("  No significant gene sets at FDR < ", fdr_cut)
      return(invisible(NULL))
    }
    
    # Build NES matrix: rows = row_id, cols = pathways
    nes_wide <- all_combined %>%
      dplyr::filter(pathway %in% sig_paths) %>%
      dplyr::select(row_id, pathway, NES) %>%
      pivot_wider(names_from = pathway, values_from = NES) %>%
      column_to_rownames("row_id")
    
    padj_wide <- all_combined %>%
      dplyr::filter(pathway %in% sig_paths) %>%
      dplyr::select(row_id, pathway, padj) %>%
      pivot_wider(names_from = pathway, values_from = padj) %>%
      column_to_rownames("row_id")
    
    # Build ordered row names: comparisons as blocks, cases within each block
    # in original case order (not clustered within blocks — dendrogram handles that)
    ordered_rows <- c()
    for (cmp in comparisons) {
      for (case_lbl in case_labels) {
        rid <- paste0(cmp$label, "  |  ", case_lbl)
        if (rid %in% rownames(nes_wide)) ordered_rows <- c(ordered_rows, rid)
      }
    }
    nes_wide  <- nes_wide[ ordered_rows,  , drop=FALSE]
    padj_wide <- padj_wide[ordered_rows,  , drop=FALSE]
    
    # Mask non-significant NES values
    for (r in rownames(nes_wide))
      for (c in colnames(nes_wide))
        if (is.na(padj_wide[r,c]) || padj_wide[r,c] >= fdr_cut)
          nes_wide[r,c] <- NA
    
    # Cluster COLUMNS across all rows
    col_mat <- nes_wide; col_mat[is.na(col_mat)] <- 0
    col_ord  <- hclust(dist(t(col_mat)))$order
    nes_wide <- nes_wide[, col_ord, drop=FALSE]
    
    # Clean column names
    colnames(nes_wide) <- str_to_title(
      gsub("_", " ", sub("^HALLMARK_", "", colnames(nes_wide)))
    )
    
    # Cluster ROWS within each comparison block separately,
    # then reassemble in block order
    clustered_rows <- c()
    for (cmp in comparisons) {
      block_rows <- ordered_rows[grepl(cmp$label, ordered_rows, fixed=TRUE)]
      block_rows <- intersect(block_rows, rownames(nes_wide))
      if (length(block_rows) < 2) {
        clustered_rows <- c(clustered_rows, block_rows)
      } else {
        bmat <- nes_wide[block_rows, , drop=FALSE]
        bmat[is.na(bmat)] <- 0
        row_clust <- hclust(dist(bmat))
        clustered_rows <- c(clustered_rows, block_rows[row_clust$order])
      }
    }
    nes_wide <- nes_wide[clustered_rows, , drop=FALSE]
    
    # Melt for ggplot
    mx <- max(abs(nes_wide), na.rm=TRUE)
    if (!is.finite(mx) || mx == 0) mx <- 1
    
    # Build comparison block boundaries for gap lines
    block_ends <- c()
    running <- 0
    for (cmp in comparisons) {
      n_in_block <- sum(grepl(cmp$label, rownames(nes_wide), fixed=TRUE))
      running <- running + n_in_block
      block_ends <- c(block_ends, running)
    }
    block_ends <- block_ends[-length(block_ends)]  # no line after last block
    
    # Row labels: strip comparison prefix, keep only case label
    row_display <- sub("^.*  \\|  ", "", rownames(nes_wide))
    # Shorten: remove "Case N: " prefix
    row_display <- gsub("^Case \\d+: ", "", row_display)
    
    # Comparison block labels (midpoint of each block)
    block_mids <- c()
    block_labels <- c()
    running <- 0
    for (cmp in comparisons) {
      n <- sum(grepl(cmp$label, rownames(nes_wide), fixed=TRUE))
      block_mids   <- c(block_mids,   running + n/2 + 0.5)
      block_labels <- c(block_labels, cmp$label)
      running <- running + n
    }
    
    hm_df <- as.data.frame(nes_wide) %>%
      rownames_to_column("row_id") %>%
      pivot_longer(-row_id, names_to="pathway", values_to="NES") %>%
      mutate(
        row_id  = factor(row_id,  levels=rev(rownames(nes_wide))),
        pathway = factor(pathway, levels=colnames(nes_wide))
      )
    
    p_hm <- ggplot(hm_df, aes(x=pathway, y=row_id, fill=NES)) +
      geom_tile(color="white", linewidth=0.3) +
      # Gap lines between comparison blocks
      geom_hline(
        yintercept = (nrow(nes_wide) - block_ends) + 0.5,
        color="grey20", linewidth=1.2
      ) +
      scale_fill_gradientn(
        colours  = c("#2980B9","#7fbbd4","white","#e87070","#C0392B"),
        values   = scales::rescale(c(-mx, -mx/3, 0, mx/3, mx)),
        limits   = c(-mx, mx),
        na.value = "grey88",
        name     = "NES"
      ) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(
        expand=c(0,0),
        labels=setNames(rev(row_display), rev(levels(hm_df$row_id)))
      ) +
      # Comparison block annotations on left
      annotate("text",
               x      = 0.2,
               y      = nrow(nes_wide) - block_mids + 1,
               label  = block_labels,
               hjust  = 1, vjust = 0.5, size = 4.5, fontface="bold",
               color  = c("#2980B9","#C0392B","#7B68EE","#2E8B57")
      ) +
      labs(
        title    = paste0("GSEA sensitivity — Classical Monocytes  (FDR < ", fdr_cut, ")"),
        subtitle = paste0("Rows = cases (8) × comparisons (4) | Columns = Hallmark gene sets\n",
                          "Blocks separated by horizontal lines | Columns clustered across all rows"),
        x=NULL, y=NULL
      ) +
      theme_minimal(base_size=14) +
      theme(
        plot.title      = element_text(face="bold", size=16, hjust=0.5),
        plot.subtitle   = element_text(size=12, hjust=0.5, color="grey40"),
        axis.text.x     = element_text(angle=45, hjust=1, size=11),
        axis.text.y     = element_text(size=11),
        legend.position = "right",
        legend.title    = element_text(size=12),
        legend.text     = element_text(size=11),
        plot.margin     = margin(t=10, r=10, b=10, l=160)
      )
    
    # Add dendrogram using ggdendro
    if (requireNamespace("ggdendro", quietly=TRUE)) {
      library(ggdendro)
      library(grid)
      library(gridExtra)
      
      # Build within-block row dendrograms stacked
      dend_grobs <- list()
      running <- 0
      total_rows <- nrow(nes_wide)
      
      for (cmp_i in seq_along(comparisons)) {
        cmp <- comparisons[[cmp_i]]
        block_rows <- rownames(nes_wide)[grepl(cmp$label, rownames(nes_wide), fixed=TRUE)]
        n <- length(block_rows)
        if (n < 2) { running <- running + n; next }
        
        bmat <- nes_wide[block_rows, , drop=FALSE]
        bmat[is.na(bmat)] <- 0
        # Rows are already in cluster order — reconstruct hclust from distance
        hc <- hclust(dist(bmat))
        dd <- dendro_data(hc, type="rectangle")
        
        # Build dendrogram plot
        p_dend <- ggplot(segment(dd)) +
          geom_segment(aes(x=x, y=y, xend=xend, yend=yend),
                       color="grey30", linewidth=0.5) +
          coord_flip() +
          scale_x_reverse() +
          theme_void() +
          theme(plot.margin=margin(0,0,0,0))
        
        dend_grobs[[cmp_i]] <- ggplotGrob(p_dend)
        running <- running + n
      }
      
      # Save with patchwork: dendrogram left + heatmap right
      if (length(dend_grobs) > 0) {
        # Use cowplot to combine if available, else just save heatmap
        if (requireNamespace("cowplot", quietly=TRUE)) {
          library(cowplot)
          dend_col <- plot_grid(plotlist=dend_grobs, ncol=1,
                                rel_heights=sapply(comparisons, function(cmp)
                                  sum(grepl(cmp$label, rownames(nes_wide), fixed=TRUE))))
          combined_full <- plot_grid(dend_col, p_hm, ncol=2,
                                     rel_widths=c(0.08, 0.92))
          ggsave(paste0("GSEA_Monocytes_combined_heatmap_FDR", fdr_cut, ".pdf"),
                 combined_full,
                 width=max(24, ncol(nes_wide)*0.4 + 10),
                 height=max(16, nrow(nes_wide)*0.4 + 4),
                 limitsize=FALSE)
          ggsave(paste0("GSEA_Monocytes_combined_heatmap_FDR", fdr_cut, ".png"),
                 combined_full,
                 width=max(24, ncol(nes_wide)*0.4 + 10),
                 height=max(16, nrow(nes_wide)*0.4 + 4),
                 dpi=120, limitsize=FALSE)
          message("  Saved combined heatmap with dendrogram (FDR < ", fdr_cut, ")")
          return(invisible(combined_full))
        }
      }
    }
    
    # Fallback: save heatmap only
    ggsave(paste0("GSEA_Monocytes_combined_heatmap_FDR", fdr_cut, ".pdf"),
           p_hm,
           width=max(24, ncol(nes_wide)*0.4 + 10),
           height=max(16, nrow(nes_wide)*0.4 + 4),
           limitsize=FALSE)
    ggsave(paste0("GSEA_Monocytes_combined_heatmap_FDR", fdr_cut, ".png"),
           p_hm,
           width=max(24, ncol(nes_wide)*0.4 + 10),
           height=max(16, nrow(nes_wide)*0.4 + 4),
           dpi=120, limitsize=FALSE)
    message("  Saved combined heatmap (FDR < ", fdr_cut, ")")
    invisible(p_hm)
  }
  
  build_combined_heatmap(FDR_LO, "0.25")
  build_combined_heatmap(FDR_HI, "0.10")
  
} else {
  message("  No CSV results found — run main loop first")
}

message("\nAll done.")


message("\nDone. Saved: GSEA_Monocytes_sensitivity_heatmap.pdf / .png")

## ── Master sensitivity heatmap (ComplexHeatmap) ─────────────────────────────
## ── 1. Load all 4 comparison CSVs ─────────────────────────────────────────────
comparison_ids <- list(
  list(id = "LentiPosNeg_LB",  label = "LentiPos vs LentiNeg (LB)"),
  list(id = "LentiPosNeg_HB",  label = "LentiPos vs LentiNeg (HB)"),
  list(id = "HBvsLB_LentiNeg", label = "HB vs LB (LentiNeg)"),
  list(id = "HBvsLB_LentiPos", label = "HB vs LB (LentiPos)")
)

case_order <- c(
  "Case 1: all cells",
  "Case 2: all cells + UMI matched",
  "Case 3: nUMI \u2265 13700",
  "Case 4: nUMI \u2265 13700 + UMI matched",
  "Case 5: nUMI \u2265 9227",
  "Case 6: nUMI \u2265 9227 + UMI matched",
  "Case 7: sample-specific SCTv2",
  "Case 8: sample-specific SCTv2 + UMI matched"
)

all_res <- bind_rows(lapply(comparison_ids, function(cmp) {
  f <- paste0("GSEA_Monocytes_", cmp$id, "_all_results.csv")
  if (!file.exists(f)) { message("Missing: ", f); return(NULL) }
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$comparison_label <- cmp$label
  df$full_id <- paste0(cmp$label, "  |  ", df$case)
  df
}))

message("Loaded ", nrow(all_res), " rows across ",
        length(unique(all_res$comparison_label)), " comparisons")

## ── 2. Filter to significant pathways ────────────────────────────────────────
sig_paths <- all_res %>%
  group_by(pathway) %>%
  summarise(min_p = min(padj, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(min_p < FDR_CUT) %>%
  pull(pathway)

message("Significant pathways at FDR < ", FDR_CUT, ": ", length(sig_paths))

if (length(sig_paths) == 0) stop("No significant pathways — try a higher FDR_CUT")

## ── 3. Build NES matrix (rows = case × comparison, cols = pathways) ──────────
mat_df <- all_res %>%
  dplyr::filter(pathway %in% sig_paths) %>%
  mutate(
    pathway_clean = str_to_title(gsub("_", " ", sub("^HALLMARK_|^MH_", "", pathway))),
    NES_masked    = ifelse(padj < FDR_CUT, NES, NA)
  ) %>%
  dplyr::select(pathway_clean, full_id, NES_masked) %>%
  pivot_wider(names_from = full_id, values_from = NES_masked) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

# Remove all-NA rows
mat_df <- mat_df[rowSums(!is.na(mat_df)) > 0, , drop = FALSE]

# Imputed matrix for clustering
mat_clust <- mat_df
mat_clust[is.na(mat_clust)] <- 0

## ── 4. Column ordering: comparison blocks × case order ───────────────────────
col_order <- c()
col_split <- c()
for (cmp in comparison_ids) {
  for (cs in case_order) {
    cid <- paste0(cmp$label, "  |  ", cs)
    if (cid %in% colnames(mat_df)) {
      col_order <- c(col_order, cid)
      col_split <- c(col_split, cmp$label)
    }
  }
}

mat_df    <- mat_df[,    col_order, drop = FALSE]
mat_clust <- mat_clust[, col_order, drop = FALSE]

# Short column labels: just the case part
col_labels <- sub("^.*  \\|  Case \\d+: ", "", colnames(mat_df))

## ── 5. Row dendrogram ─────────────────────────────────────────────────────────
row_dend <- hclust(dist(mat_clust, method = "euclidean"), method = "ward.D2")

## ── 6. Colour scale ───────────────────────────────────────────────────────────
max_nes <- max(abs(mat_df), na.rm = TRUE)
col_fun <- colorRamp2(
  c(-max_nes, -max_nes/2, 0, max_nes/2, max_nes),
  c("#2980B9", "#7fbbd4", "white", "#e87070", "#C0392B")
)

## ── 7. Draw heatmap ───────────────────────────────────────────────────────────
ht <- Heatmap(
  matrix          = mat_df,
  col             = col_fun,
  na_col          = "grey92",
  name            = "NES",
  
  # Row clustering
  cluster_rows    = row_dend,
  show_row_dend   = TRUE,
  row_dend_width  = unit(4, "cm"),
  row_names_side  = "left",
  row_names_gp    = gpar(fontsize = 11),
  row_names_max_width = unit(20, "cm"),
  
  # Column ordering (no clustering — blocks fixed)
  cluster_columns  = FALSE,
  column_split     = factor(col_split, levels = sapply(comparison_ids, `[[`, "label")),
  column_gap       = unit(8, "mm"),
  column_labels    = col_labels,
  column_names_gp  = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_max_height = unit(8, "cm"),
  column_title_gp  = gpar(fontsize = 13, fontface = "bold"),
  
  # Cell borders
  rect_gp = gpar(col = "white", lwd = 1.2),
  
  # Legend
  heatmap_legend_param = list(
    title_gp   = gpar(fontsize = 12, fontface = "bold"),
    labels_gp  = gpar(fontsize = 11),
    legend_height = unit(6, "cm"),
    grid_width    = unit(0.8, "cm")
  )
)

## ── 8. Save ───────────────────────────────────────────────────────────────────
out_stem <- paste0("GSEA_Master_Sensitivity_Mono_FDR",
                   gsub("\\.", "_", as.character(FDR_CUT)))

pdf(paste0(out_stem, ".pdf"),
    width  = max(20, ncol(mat_df) * 0.55 + 10),
    height = max(12, nrow(mat_df) * 0.35 + 4))

draw(ht,
     column_title    = paste0("Classical Monocytes — Global GSEA Sensitivity: Hallmark Pathways (FDR < ", FDR_CUT, ")"),
     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
     padding         = unit(c(3, 3, 3, 6), "cm"),
     merge_legend    = TRUE)

dev.off()

message("Saved: ", out_stem, ".pdf")

