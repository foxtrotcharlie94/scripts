## ============================================================
##  GSEA sensitivity analysis
##  4 comparisons x 4 UMI-filtering cases x 2 FDR thresholds
##  Output: 4x2 grid of Hallmark heatmaps
## ============================================================

library(data.table)
library(tidyverse)
library(edgeR)
library(fgsea)
library(msigdbr)
library(patchwork)

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
N_PERM       <- 1000

## ── Load Hallmark gene sets ───────────────────────────────────────────────────
message("Loading Hallmark gene sets...")
gs_hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
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
  ranked <- setNames(sign(deg$logFC) * (-log10(deg$PValue)), deg$gene)
  ranked <- ranked[!duplicated(names(ranked)) & is.finite(ranked)]
  fgsea(pathways = gs_hallmark, stats = ranked,
        minSize = 15, maxSize = 500, nPermSimple = N_PERM) %>%
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

message("\nDone. Saved: GSEA_Monocytes_sensitivity_heatmap.pdf / .png")

## =============================================================================
##  GSEA Master Sensitivity Analysis: Integrated Clustered Heatmap
##  Output: A 4-Comparison Faceted Heatmap with Row Dendrograms
## =============================================================================

library(data.table)
library(tidyverse)
library(edgeR)
library(fgsea)
library(msigdbr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# Ensure results are loaded from previous steps (if running in parts)
# all_combined_res <- bind_rows(lapply(comparisons, function(c) {
#   df <- read.csv(paste0("GSEA_Monocytes_", c$id, "_all_results.csv"))
#   df$comparison_label <- c$label
#   df$full_id <- paste(df$comparison_label, df$case, sep = " | ")
#   return(df)
# }))

message("\n=== Building Refined Master Heatmap with FDR Masking ===")

## ── 1. Advanced Filtering ───────────────────────────────────────────────────
# Keep pathways that are significant (FDR < 0.05) in AT LEAST ONE case/comparison
sig_pathways <- all_combined_res %>%
  group_by(pathway) %>%
  summarise(min_p = min(padj, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(min_p < 0.05) %>%
  pull(pathway)

## ── 2. Build Matrix with FDR Masking ──────────────────────────────────────────
# We replace NES with NA if the FDR threshold isn't met to trigger the grey-out
master_mat_raw <- all_combined_res %>%
  dplyr::filter(pathway %in% sig_pathways) %>%
  mutate(
    # Clean up pathway names for the plot
    pathway_clean = str_to_title(gsub("_", " ", sub("^HALLMARK_", "", pathway))),
    # Masking: Only keep NES if FDR < 0.05
    NES_masked = ifelse(padj < 0.05, NES, NA)
  ) %>%
  select(pathway_clean, full_id, NES_masked) %>%
  pivot_wider(names_from = full_id, values_from = NES_masked) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

## ── 3. Build Matrix & Prepare for Clustering ────────────────────────────────
# Step A: Create the visual matrix (with NAs for the grey-out)
master_mat_visual <- all_combined_res %>%
  dplyr::filter(pathway %in% sig_pathways) %>%
  mutate(
    pathway_clean = str_to_title(gsub("_", " ", sub("^HALLMARK_", "", pathway))),
    NES_masked = ifelse(padj < 0.05, NES, NA)
  ) %>%
  select(pathway_clean, full_id, NES_masked) %>%
  pivot_wider(names_from = full_id, values_from = NES_masked) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

# Step B: Filter out rows that are 100% NA (no signal in any case)
keep_rows <- rowSums(!is.na(master_mat_visual)) > 0
master_mat_visual <- master_mat_visual[keep_rows, , drop = FALSE]

# Step C: Create the clustering matrix (Impute NAs to 0)
# This is ONLY used to calculate the dendrogram
master_mat_cluster <- master_mat_visual
master_mat_cluster[is.na(master_mat_cluster)] <- 0

## ── 4. Generate the Heatmap (Maximum Space & Large Font Edition) ────────────
# Pre-calculate the dendrogram using the imputed matrix
row_dend = hclust(dist(master_mat_cluster, method = "euclidean"), method = "ward.D2")

# Shorten bottom labels to the absolute minimum to prevent overlap
clean_column_labels <- sub(".*Case ", "Case ", colnames(master_mat_visual))

final_plot <- Heatmap(
  matrix = master_mat_visual,
  cluster_rows = row_dend,
  
  # --- Color & Scale ---
  name = "NES",
  col = col_fun,
  na_col = "grey95",
  
  # --- Heatmap Dimensions ---
  # Increase the body width significantly
  width = unit(70, "cm"), 
  
  # --- Top Titles (Comparison Facets) ---
  cluster_columns = FALSE,
  column_split = col_meta$comparison, 
  column_gap = unit(20, "mm"), # DOUBLED gap to prevent title overlap
  
  # Styling the comparison names (HB vs LB, etc.)
  column_title_gp = gpar(fontsize = 34, fontface = "bold"), 
  
  # --- Bottom Labels (Case Numbers) ---
  column_labels = clean_column_labels,
  column_names_gp = gpar(fontsize = 27, fontface = "plain"), 
  column_names_rot = 45,
  column_names_max_height = unit(12, "cm"),
  
  # --- Side Labels (Pathway Names) ---
  show_row_dend = TRUE,
  row_dend_width = unit(6, "cm"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 28, fontface = "bold"), # Much larger pathways
  row_names_max_width = unit(30, "cm"), 
  
  # --- Legend Configuration ---
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 22, fontface = "bold"),
    labels_gp = gpar(fontsize = 20),
    legend_height = unit(12, "cm"),
    grid_width = unit(1.5, "cm")
  ),
  
  # --- Aesthetics ---
  rect_gp = gpar(col = "white", lwd = 2) # Thicker grid lines for large scale
)

# 5. Save to a Massive Canvas
# Using 50x35 inches gives R an enormous amount of physical space
pdf("GSEA_Master_Sensitivity_Ultra_HD.pdf", width = 50, height = 35)

draw(final_plot, 
     merge_legend = TRUE, 
     # Main Overarching Title
     column_title = "Classical Monocytes - Global GSEA Sensitivity: Hallmark Pathways (FDR < 0.05)",
     column_title_gp = gpar(fontsize = 36, fontface = "bold"),
     # Massive padding to utilize the 50-inch width
     # Top, Right, Bottom, Left
     padding = unit(c(5, 5, 15, 25), "cm")) 

dev.off()

