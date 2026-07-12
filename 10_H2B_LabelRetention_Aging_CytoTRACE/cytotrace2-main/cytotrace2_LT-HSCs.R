
library(CytoTRACE2)

## ============================================================
## CytoTRACE — differentiation state scoring
##
## Runs on the sub5000 Seurat object (or swap in the fullUMI
## one by changing the RDS path).
##
## CytoTRACE score: 1 = most stem-like, 0 = most differentiated
##
## Outputs:
##   CytoTRACE_UMAP_score.png
##   CytoTRACE_UMAP_gcs.png        (gene count score)
##   CytoTRACE_violin_group.png
##   CytoTRACE_violin_sample.png
##   CytoTRACE_vs_pseudotime.png   (correlation with Slingshot)
##   CytoTRACE_stats.txt
## ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(lme4)
library(patchwork)

set.seed(42)

OUTDIR <- "C:/Users/fc809/Downloads/Sub5000_output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)



## ── 2. LOAD SEURAT OBJECT ─────────────────────────────────────────────────────

cat("Loading Seurat object...\n")
so <- readRDS(file.path(OUTDIR, "LT_HSC_sub5000_seurat.rds"))
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(so), nrow(so)))
cat("Groups:\n"); print(table(so$group))

## ── 3. RUN CytoTRACE2 ─────────────────────────────────────────────────────────
# Input: log-normalised counts matrix (genes x cells)
# species: "mouse" for GRCm39

cat("Running CytoTRACE2...\n")

expr_mat <- as.matrix(GetAssayData(so, layer = "data"))  # log-normalised

cytotrace_result <- cytotrace2(
  input   = GetAssayData(so, layer = "counts"),
  species = "mouse",
  ncores  = 1
)


## ── 4. ADD TO SEURAT METADATA ─────────────────────────────────────────────────
# cytotrace_result is a dataframe with cells as rownames

so$cytotrace_score    <- cytotrace_result[colnames(so), "CytoTRACE2_Score"]
so$cytotrace_potency  <- cytotrace_result[colnames(so), "CytoTRACE2_Potency"]
so$cytotrace_relative <- cytotrace_result[colnames(so), "CytoTRACE2_Relative"]
cat("Non-NA scores:", sum(!is.na(so$cytotrace_score)), "\n")

## ── 5. PLOTS ──────────────────────────────────────────────────────────────────

group_cols <- c(LB_LentiNeg="#90CAF9", LB_LentiPos="#1A237E",
                HB_LentiNeg="#EF9F27", HB_LentiPos="#B71C1C")

umap_df              <- as.data.frame(Embeddings(so, "umap"))
colnames(umap_df)    <- c("UMAP1","UMAP2")
umap_df$ct_score     <- so$cytotrace_score
umap_df$ct_relative  <- so$cytotrace_relative
umap_df$group        <- so$group
umap_df$sample       <- so$orig.ident
umap_df$division     <- so$division_burden

# Add pseudotime safely if it exists
umap_df$pseudotime <- tryCatch(so$pseudotime_ss, error = function(e) NA_real_)

# UMAP coloured by CytoTRACE2 score
p_ct_umap <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=ct_score)) +
  geom_point(size=0.5, alpha=0.7) +
  scale_color_viridis_c(option="magma", direction=-1,
                        name="CytoTRACE2\nscore\n(1=stem)") +
  labs(title="CytoTRACE2 score") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0.5))

# UMAP coloured by relative order
p_rel_umap <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=ct_relative)) +
  geom_point(size=0.5, alpha=0.7) +
  scale_color_viridis_c(option="cividis", direction=-1,
                        name="Relative\norder\n(1=stem)") +
  labs(title="CytoTRACE2 relative order") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0.5))

ggsave(file.path(OUTDIR, "CytoTRACE_UMAP_score.png"),
       p_ct_umap, width=6, height=5, dpi=150)
ggsave(file.path(OUTDIR, "CytoTRACE_UMAP_relative.png"),
       p_rel_umap, width=6, height=5, dpi=150)
message("Saved: CytoTRACE UMAPs")

# Violin by group
vln_df       <- umap_df
vln_df$group <- factor(vln_df$group,
                       levels=c("LB_LentiNeg","LB_LentiPos",
                                "HB_LentiNeg","HB_LentiPos"))

p_vln_grp <- ggplot(vln_df, aes(x=group, y=ct_score, fill=group)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.8) +
  geom_boxplot(width=0.12, outlier.size=0.3, fill="white", alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=23, size=2.5,
               fill="white", color="black") +
  scale_fill_manual(values=group_cols, guide="none") +
  labs(x=NULL, y="CytoTRACE2 score (1=stem)",
       title="Differentiation state by group") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title=element_text(hjust=0.5))

# Violin by sample
p_vln_samp <- ggplot(vln_df, aes(x=sample, y=ct_score, fill=division)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.8) +
  geom_boxplot(width=0.12, outlier.size=0.3, fill="white", alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=23, size=2.5,
               fill="white", color="black") +
  scale_fill_manual(values=c(HB="#B71C1C", LB="#1A237E"), name="Burden") +
  labs(x=NULL, y="CytoTRACE2 score (1=stem)",
       title="Differentiation state by sample") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title=element_text(hjust=0.5))

ggsave(file.path(OUTDIR, "CytoTRACE_violin_group.png"),
       p_vln_grp, width=6, height=5, dpi=150)
ggsave(file.path(OUTDIR, "CytoTRACE_violin_sample.png"),
       p_vln_samp, width=8, height=5, dpi=150)
message("Saved: CytoTRACE violins")

# Correlation with Slingshot pseudotime (if available)
pt_available <- !all(is.na(umap_df$pseudotime)) &&
  sum(!is.na(umap_df$pseudotime)) > 100

if (pt_available) {
  cor_df <- umap_df[!is.na(umap_df$pseudotime) & !is.na(umap_df$ct_score), ]
  r_val  <- round(cor(cor_df$pseudotime, cor_df$ct_score,
                      method="spearman"), 3)
  p_corr <- ggplot(cor_df, aes(x=pseudotime, y=ct_score, color=group)) +
    geom_point(size=0.4, alpha=0.5) +
    geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.8) +
    scale_color_manual(values=group_cols, name="Group") +
    annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5, size=3.5,
             label=paste0("Spearman r = ", r_val)) +
    labs(x="Slingshot pseudotime", y="CytoTRACE2 score (1=stem)",
         title="CytoTRACE2 vs Slingshot pseudotime") +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(hjust=0.5))
  ggsave(file.path(OUTDIR, "CytoTRACE_vs_pseudotime.png"),
         p_corr, width=7, height=5, dpi=150)
  message(sprintf("Saved: CytoTRACE_vs_pseudotime.png  (Spearman r = %s)", r_val))
} else {
  cat("Pseudotime not in metadata — skipping correlation plot.\n")
  cat("Add so$pseudotime_ss from a Slingshot run and rerun this block.\n")
}

## ── 6. STATS ──────────────────────────────────────────────────────────────────

lm_df <- data.frame(
  ct_score = so$cytotrace_score,
  HB       = as.integer(so$division_burden == "HB"),
  LentiPos = as.integer(so$lenti_status    == "LentiPos"),
  sample   = so$orig.ident
)
lm_df <- lm_df[!is.na(lm_df$ct_score), ]
cat("Rows in lm_df:", nrow(lm_df), "\n")

# LentiPos — mixed model (within-animal)
lme_lenti <- lmer(ct_score ~ LentiPos + (1 | sample), data = lm_df)
cat("\n=== LentiPos mixed model (CytoTRACE2 score) ===\n")
print(summary(lme_lenti))

# HB/LB — animal-level t-test
animal_means <- aggregate(mean_ct ~ sample,
                          data = data.frame(mean_ct = lm_df$ct_score,
                                            sample  = lm_df$sample),
                          FUN = mean)
animal_means$HB <- factor(ifelse(grepl("HB", animal_means$sample), 1, 0),
                          levels = c(0, 1))

t_hb <- t.test(mean_ct ~ HB, data = animal_means)
cat("\n=== HB t-test (animal-level CytoTRACE2 means) ===\n")
print(t_hb)
cat("\nAnimal means:\n"); print(animal_means)

# Save stats
sink(file.path(OUTDIR, "CytoTRACE_stats.txt"))
cat("=== CytoTRACE2 — sub5000 Seurat object ===\n\n")
cat("Score: 1 = most stem-like, 0 = most differentiated\n\n")
cat("Potency breakdown:\n")
print(table(so$cytotrace_potency))
cat("\nLentiPos mixed model: ct_score ~ LentiPos + (1|sample)\n\n")
print(summary(lme_lenti))
cat("\nHB t-test (animal-level CytoTRACE2 means):\n")
print(t_hb)
cat("\nAnimal means:\n"); print(as.data.frame(animal_means))
if (pt_available) cat(sprintf("\nSpearman r vs Slingshot pseudotime: %s\n", r_val))
sink()
message("Saved: CytoTRACE_stats.txt")

## ── 7. SAVE UPDATED SEURAT OBJECT ─────────────────────────────────────────────

saveRDS(so, file.path(OUTDIR, "LT_HSC_sub5000_seurat.rds"))
message("Seurat object updated with CytoTRACE2 scores and resaved.")

message("\n=== Done. Outputs in: ", OUTDIR, " ===")

## ============================================================
## Combined publication figure: CytoTRACE2 + Slingshot
##
## Assumes so already has so$cytotrace_score from CytoTRACE2_analysis.R
## Re-runs pooled Slingshot on clusters 1-5 to get pseudotime,
## then builds a 4-panel publication figure:
##
##   Panel A: UMAP coloured by CytoTRACE2 score
##   Panel B: UMAP coloured by Slingshot pseudotime
##   Panel C: Violin — CytoTRACE2 score by group
##   Panel D: Violin — Slingshot pseudotime by group
##
## Output: Combined_figure.png + Combined_figure_violins_only.png
## ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(slingshot)
library(SingleCellExperiment)

set.seed(42)

OUTDIR <- "C:/Users/fc809/Downloads/Sub5000_output"

group_cols <- c(LB_LentiNeg="#90CAF9", LB_LentiPos="#1A237E",
                HB_LentiNeg="#EF9F27", HB_LentiPos="#B71C1C")

## ── 1. LOAD SEURAT OBJECT ─────────────────────────────────────────────────────

cat("Loading Seurat object...\n")
so <- readRDS(file.path(OUTDIR, "LT_HSC_sub5000_seurat.rds"))
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(so), nrow(so)))

if (is.null(so$cytotrace_score) || all(is.na(so$cytotrace_score))) {
  stop("CytoTRACE2 scores not found. Run CytoTRACE2_analysis.R first.")
}
cat("CytoTRACE2 scores found:", sum(!is.na(so$cytotrace_score)), "cells\n")

## ── 2. POOLED SLINGSHOT — clusters 1-5, all samples ──────────────────────────

cat("Running pooled Slingshot on clusters 1-5...\n")

so_c15 <- subset(so, subset = seurat_clusters %in% c(1,2,3,4,5))
cat(sprintf("Clusters 1-5: %d cells\n", ncol(so_c15)))

sce <- as.SingleCellExperiment(so_c15)
reducedDim(sce, "UMAP") <- Embeddings(so_c15, "umap")
sce$division_burden     <- so_c15$division_burden
sce$lenti_status        <- so_c15$lenti_status
sce$group               <- so_c15$group
sce$sling_cluster       <- factor(rep("all", ncol(sce)))

# Root: closest cell to LB/LentiNeg centroid
umap_coords <- reducedDim(sce, "UMAP")
lb_lneg_idx <- which(sce$division_burden == "LB" & sce$lenti_status == "LentiNeg")
centroid    <- colMeans(umap_coords[lb_lneg_idx, , drop=FALSE])
dists       <- sqrt(rowSums(sweep(umap_coords, 2, centroid)^2))

sce <- slingshot(sce, clusterLabels="sling_cluster",
                 reducedDim="UMAP", start.clus="all", approx_points=150)

pt_raw  <- slingPseudotime(sce)[, 1]
mean_lb <- mean(pt_raw[lb_lneg_idx], na.rm=TRUE)
mean_hb <- mean(pt_raw[which(sce$division_burden == "HB")], na.rm=TRUE)
if (mean_lb > mean_hb) pt_raw <- max(pt_raw, na.rm=TRUE) - pt_raw

# Write pseudotime back to so (only for cells in clusters 1-5)
so$pseudotime_ss <- NA_real_
so$pseudotime_ss[colnames(so_c15)] <- pt_raw
cat(sprintf("Pseudotime range: %.2f — %.2f\n",
            min(pt_raw, na.rm=TRUE), max(pt_raw, na.rm=TRUE)))

curve_df           <- as.data.frame(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ])
colnames(curve_df) <- c("UMAP1","UMAP2")

## ── 3. BUILD PLOT DATA ────────────────────────────────────────────────────────

umap_df           <- as.data.frame(Embeddings(so, "umap"))
colnames(umap_df) <- c("UMAP1","UMAP2")
umap_df$ct_score  <- so$cytotrace_score
umap_df$pseudotime <- so$pseudotime_ss
umap_df$group     <- factor(so$group,
                            levels=c("LB_LentiNeg","LB_LentiPos",
                                     "HB_LentiNeg","HB_LentiPos"))
umap_df$division  <- so$division_burden

## ── 4. PANEL A — UMAP: CytoTRACE2 score ─────────────────────────────────────

pA <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=ct_score)) +
  geom_point(size=0.4, alpha=0.7) +
  scale_color_viridis_c(option="magma", direction=-1,
                        name="CytoTRACE2\nscore") +
  labs(title="A  Differentiation state (CytoTRACE2)") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0, face="bold"),
        legend.title=element_text(size=9))

## ── 5. PANEL B — UMAP: Slingshot pseudotime ──────────────────────────────────

umap_pt <- umap_df[!is.na(umap_df$pseudotime), ]

pB <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=pseudotime), size=0.4, alpha=0.7) +
  geom_path(data=curve_df, aes(x=UMAP1, y=UMAP2),
            color="black", linewidth=0.7, inherit.aes=FALSE) +
  scale_color_viridis_c(option="plasma", name="Pseudotime",
                        na.value="grey85") +
  labs(title="B  Trajectory (Slingshot, clusters 1–5)") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0, face="bold"),
        legend.title=element_text(size=9))

## ── 6. PANEL C — Violin: CytoTRACE2 by group ─────────────────────────────────

pC <- ggplot(umap_df, aes(x=group, y=ct_score, fill=group)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.85) +
  geom_boxplot(width=0.1, outlier.size=0.3, fill="white", alpha=0.8) +
  stat_summary(fun=mean, geom="point", shape=23, size=2.5,
               fill="white", color="black") +
  scale_fill_manual(values=group_cols, guide="none") +
  labs(x=NULL, y="CytoTRACE2 score (1=stem)",
       title="C  Differentiation state by group") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title=element_text(hjust=0, face="bold"))

## ── 7. PANEL D — Violin: Slingshot pseudotime by group ───────────────────────

vln_pt       <- umap_pt
vln_pt$group <- factor(vln_pt$group,
                       levels=c("LB_LentiNeg","LB_LentiPos",
                                "HB_LentiNeg","HB_LentiPos"))
grp_counts   <- table(vln_pt$group)
vln_pt       <- vln_pt[vln_pt$group %in% names(grp_counts[grp_counts >= 3]), ]
vln_pt$group <- droplevels(vln_pt$group)

pD <- ggplot(vln_pt, aes(x=group, y=pseudotime, fill=group)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.85) +
  geom_boxplot(width=0.1, outlier.size=0.3, fill="white", alpha=0.8) +
  stat_summary(fun=mean, geom="point", shape=23, size=2.5,
               fill="white", color="black") +
  scale_fill_manual(values=group_cols, guide="none", drop=TRUE) +
  labs(x=NULL, y="Pseudotime (early → late)",
       title="D  Pseudotime by group") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title=element_text(hjust=0, face="bold"))

## ── 8. COMBINE AND SAVE ───────────────────────────────────────────────────────

# 4-panel figure
p_combined <- (pA | pB) / (pC | pD)
ggsave(file.path(OUTDIR, "Combined_figure.png"),
       p_combined, width=12, height=10, dpi=200)
message("Saved: Combined_figure.png")

# Violins only (cleaner for publication)
p_violins <- pC | pD
ggsave(file.path(OUTDIR, "Combined_figure_violins_only.png"),
       p_violins, width=10, height=5, dpi=200)
message("Saved: Combined_figure_violins_only.png")

# Save updated Seurat object with pseudotime
saveRDS(so, file.path(OUTDIR, "LT_HSC_sub5000_seurat.rds"))
message("Seurat object updated with pseudotime_ss and resaved.")

message("\n=== Done. Outputs in: ", OUTDIR, " ===")

## ============================================================
## CytoTRACE2 — effect size contextualization
##
## Computes multiple ways of expressing the HB/LB and LentiPos
## effect sizes to find the most interpretable framing.
##
## Requires: so with cytotrace_score and cytotrace_relative
## install.packages("effsize") if needed
## ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(lme4)
library(effsize)

set.seed(42)

OUTDIR <- "C:/Users/fc809/Downloads/Sub5000_output"
so <- readRDS(file.path(OUTDIR, "LT_HSC_sub5000_seurat.rds"))

## ── Build working dataframe ───────────────────────────────────────────────────

df <- data.frame(
  ct_score    = so$cytotrace_score,
  ct_relative = so$cytotrace_relative,
  HB          = so$division_burden == "HB",
  LentiPos    = so$lenti_status == "LentiPos",
  sample      = as.character(so$orig.ident),
  group       = so$group,
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$ct_score), ]
cat(sprintf("Cells in analysis: %d\n", nrow(df)))

## ── 1. RAW DIFFERENCE ─────────────────────────────────────────────────────────

hb_mean  <- mean(df$ct_score[df$HB])
lb_mean  <- mean(df$ct_score[!df$HB])
raw_diff <- hb_mean - lb_mean   # negative = HB more differentiated

cat("\n══════════════════════════════════════════════\n")
cat("  1. Raw score difference\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  LB mean:  %.5f\n", lb_mean))
cat(sprintf("  HB mean:  %.5f\n", hb_mean))
cat(sprintf("  HB - LB:  %.5f\n", raw_diff))

## ── 2. % OF OBSERVED RANGE ────────────────────────────────────────────────────

obs_range <- diff(range(df$ct_score))
pct_range <- abs(raw_diff) / obs_range * 100

cat("\n══════════════════════════════════════════════\n")
cat("  2. % of observed score range\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  Observed range: %.3f — %.3f  (width = %.3f)\n",
            min(df$ct_score), max(df$ct_score), obs_range))
cat(sprintf("  Effect = %.1f%% of observed range\n", pct_range))

## ── 3. % OF OLIGOPOTENT BIN WIDTH ─────────────────────────────────────────────
# CytoTRACE2 bins: Differentiated 0-0.1, Unipotent 0.1-0.2,
#                  Oligopotent 0.2-0.4, Multipotent 0.4-0.6,
#                  Pluripotent 0.6-0.8, Totipotent 0.8-1.0

oligo_width <- 0.2
pct_oligo   <- abs(raw_diff) / oligo_width * 100

cat("\n══════════════════════════════════════════════\n")
cat("  3. % of Oligopotent bin width (0.2–0.4)\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  Effect = %.1f%% of Oligopotent bin\n", pct_oligo))

## ── 4. COHEN'S D (cell-level) ─────────────────────────────────────────────────

cd_hb <- cohen.d(ct_score ~ factor(ifelse(HB, "HB", "LB")), data=df)

cat("\n══════════════════════════════════════════════\n")
cat("  4. Cohen's d — HB vs LB (cell-level)\n")
cat("══════════════════════════════════════════════\n")
print(cd_hb)

## ── 5. COHEN'S D (animal-level) ───────────────────────────────────────────────

animal_means <- aggregate(ct_score ~ sample, data=df, FUN=mean)
animal_means$HB <- ifelse(grepl("HB", animal_means$sample), "HB", "LB")
cd_hb_animal <- cohen.d(ct_score ~ factor(HB), data=animal_means)

cat("\n══════════════════════════════════════════════\n")
cat("  5. Cohen's d — HB vs LB (animal-level means)\n")
cat("══════════════════════════════════════════════\n")
print(cd_hb_animal)

## ── 6. CytoTRACE2_RELATIVE — rank-normalized score ───────────────────────────

hb_rel   <- mean(df$ct_relative[df$HB],  na.rm=TRUE)
lb_rel   <- mean(df$ct_relative[!df$HB], na.rm=TRUE)
rel_diff <- hb_rel - lb_rel

cat("\n══════════════════════════════════════════════\n")
cat("  6. CytoTRACE2_Relative (rank-normalized 0-1)\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  LB mean relative: %.4f\n", lb_rel))
cat(sprintf("  HB mean relative: %.4f\n", hb_rel))
cat(sprintf("  HB - LB:          %.4f  (= %.1f percentile points)\n",
            rel_diff, abs(rel_diff)*100))

# Animal-level t-test on relative score
animal_rel     <- aggregate(ct_relative ~ sample, data=df, FUN=mean)
animal_rel$HB  <- factor(ifelse(grepl("HB", animal_rel$sample), 1, 0),
                         levels=c(0,1))
t_rel <- t.test(ct_relative ~ HB, data=animal_rel)
cat("\nAnimal-level t-test on CytoTRACE2_Relative:\n")
print(t_rel)

# Mixed model on relative score — LentiPos effect in percentile points
lme_rel  <- lmer(ct_relative ~ LentiPos + (1|sample), data=df)
coef_rel <- coef(summary(lme_rel))
cat("\nLentiPos mixed model on CytoTRACE2_Relative:\n")
print(coef_rel)
lenti_row  <- rownames(coef_rel)[grepl("LentiPos", rownames(coef_rel))]
lenti_beta <- coef_rel[lenti_row, "Estimate"]
lenti_t    <- coef_rel[lenti_row, "t value"]
cat(sprintf("\n  LentiPos beta = %.5f, t = %.3f\n", lenti_beta, lenti_t))
cat(sprintf("  LentiPos effect = %.1f percentile points\n",
            abs(lenti_beta)*100))

## ── 7. EFFECT RELATIVE TO BETWEEN-CATEGORY DISTANCE ─────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  7. Effect vs between-category distance\n")
cat("══════════════════════════════════════════════\n")

cat_means <- tapply(df$ct_score,
                    as.character(so$cytotrace_potency[!is.na(so$cytotrace_score)]),
                    mean)
cat("Mean score per potency category:\n")
print(round(sort(cat_means), 4))

unipotent_mean      <- cat_means["Unipotent"]
oligo_mean          <- cat_means["Oligopotent"]
multi_mean          <- cat_means["Multipotent"]
unipotent_oligo_gap <- abs(oligo_mean - unipotent_mean)
oligo_multi_gap     <- abs(multi_mean - oligo_mean)

cat(sprintf("\n  Unipotent → Oligopotent gap:   %.4f\n", unipotent_oligo_gap))
cat(sprintf("  Oligopotent → Multipotent gap: %.4f\n", oligo_multi_gap))
cat(sprintf("\n  HB-LB as %% of Unipotent→Oligo gap: %.1f%%\n",
            abs(raw_diff)/unipotent_oligo_gap*100))
cat(sprintf("  HB-LB as %% of Oligo→Multi gap:     %.1f%%\n",
            abs(raw_diff)/oligo_multi_gap*100))

## ── 8. SUMMARY TABLE ──────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  SUMMARY — HB vs LB effect size\n")
cat("══════════════════════════════════════════════\n")

summary_tbl <- data.frame(
  Metric = c(
    "Raw score difference (HB-LB)",
    "% of observed score range",
    "% of Oligopotent bin width",
    "Cohen's d (cell-level)",
    "Cohen's d (animal-level)",
    "Relative score difference (percentile pts)",
    "% of Unipotent->Oligo category gap",
    "% of Oligo->Multipotent category gap"
  ),
  Value = c(
    round(raw_diff, 5),
    round(pct_range, 2),
    round(pct_oligo, 2),
    round(cd_hb$estimate, 4),
    round(cd_hb_animal$estimate, 4),
    round(abs(rel_diff)*100, 2),
    round(abs(raw_diff)/unipotent_oligo_gap*100, 1),
    round(abs(raw_diff)/oligo_multi_gap*100, 1)
  ),
  stringsAsFactors = FALSE
)
print(summary_tbl, row.names=FALSE)

write.csv(summary_tbl,
          file.path(OUTDIR, "CytoTRACE_effect_size_summary.csv"),
          row.names=FALSE)
message("Saved: CytoTRACE_effect_size_summary.csv")

# One mean per animal x group combination
compartment_means <- aggregate(ct_score ~ sample + group, data=df, FUN=mean)
compartment_means$HB <- ifelse(grepl("HB", compartment_means$sample), "HB", "LB")

cat("Compartment means:\n")
print(compartment_means)

# Cohen's d on HB vs LB at compartment level
cd_compartment <- cohen.d(ct_score ~ factor(HB), data=compartment_means)
cat("\nCohen's d — HB vs LB (compartment-level means):\n")
print(cd_compartment)

# Also compute for LentiPos vs LentiNeg at compartment level
compartment_means$LentiPos <- ifelse(grepl("LentiPos", compartment_means$group),
                                     "LentiPos", "LentiNeg")
cd_lenti_compartment <- cohen.d(ct_score ~ factor(LentiPos),
                                data=compartment_means)
cat("\nCohen's d — LentiPos vs LentiNeg (compartment-level means):\n")
print(cd_lenti_compartment)

# Cohen's d per potency category — HB vs LB
cat("══════════════════════════════════════════════\n")
cat("  Cohen's d by potency category — HB vs LB\n")
cat("══════════════════════════════════════════════\n")

df$potency <- as.character(so$cytotrace_potency[!is.na(so$cytotrace_score)])

categories <- names(table(df$potency)[table(df$potency) >= 10])

for (cat_name in categories) {
  df_cat <- df[df$potency == cat_name, ]
  n_hb   <- sum(df_cat$HB)
  n_lb   <- sum(!df_cat$HB)
  
  if (n_hb < 5 | n_lb < 5) {
    cat(sprintf("\n  %-15s — skipped (HB n=%d, LB n=%d)\n",
                cat_name, n_hb, n_lb))
    next
  }
  
  cd <- tryCatch(
    cohen.d(ct_score ~ factor(ifelse(HB, "HB", "LB")), data=df_cat),
    error = function(e) NULL
  )
  
  if (!is.null(cd)) {
    cat(sprintf("\n  %-15s  n=%d (HB=%d, LB=%d)\n",
                cat_name, nrow(df_cat), n_hb, n_lb))
    cat(sprintf("    d = %.4f (%s)\n", cd$estimate, cd$magnitude))
    cat(sprintf("    95%% CI: [%.3f, %.3f]\n",
                cd$conf.int[1], cd$conf.int[2]))
  }
}

# Cohen's d per potency category — HB vs LB
cat("══════════════════════════════════════════════\n")
cat("  Cohen's d by potency category — HB vs LB\n")
cat("══════════════════════════════════════════════\n")

df$potency <- as.character(so$cytotrace_potency[!is.na(so$cytotrace_score)])

categories <- names(table(df$potency)[table(df$potency) >= 10])

for (cat_name in categories) {
  df_cat <- df[df$potency == cat_name, ]
  n_hb   <- sum(df_cat$HB)
  n_lb   <- sum(!df_cat$HB)
  
  if (n_hb < 5 | n_lb < 5) {
    cat(sprintf("\n  %-15s — skipped (HB n=%d, LB n=%d)\n",
                cat_name, n_hb, n_lb))
    next
  }
  
  cd <- tryCatch(
    cohen.d(ct_score ~ factor(ifelse(HB, "HB", "LB")), data=df_cat),
    error = function(e) NULL
  )
  
  if (!is.null(cd)) {
    cat(sprintf("\n  %-15s  n=%d (HB=%d, LB=%d)\n",
                cat_name, nrow(df_cat), n_hb, n_lb))
    cat(sprintf("    d = %.4f (%s)\n", cd$estimate, cd$magnitude))
    cat(sprintf("    95%% CI: [%.3f, %.3f]\n",
                cd$conf.int[1], cd$conf.int[2]))
  }
}

# Also compartment-level per category
cat("\n══════════════════════════════════════════════\n")
cat("  Cohen's d by potency category — compartment level\n")
cat("══════════════════════════════════════════════\n")

for (cat_name in categories) {
  df_cat <- df[df$potency == cat_name, ]
  n_hb   <- sum(df_cat$HB)
  n_lb   <- sum(!df_cat$HB)
  if (n_hb < 5 | n_lb < 5) next
  
  comp <- aggregate(ct_score ~ sample + group, data=df_cat, FUN=mean)
  comp$HB <- ifelse(grepl("HB", comp$sample), "HB", "LB")
  
  if (length(unique(comp$HB)) < 2) next
  
  cd <- tryCatch(
    cohen.d(ct_score ~ factor(HB), data=comp),
    error = function(e) NULL
  )
  
  if (!is.null(cd)) {
    cat(sprintf("\n  %-15s  n_compartments=%d\n", cat_name, nrow(comp)))
    cat(sprintf("    d = %.4f (%s)\n", cd$estimate, cd$magnitude))
    cat(sprintf("    95%% CI: [%.3f, %.3f]\n",
                cd$conf.int[1], cd$conf.int[2]))
  }
}

# ── Animal level (all cells) ──────────────────────────────────────────────────
animal_lenti_all <- aggregate(ct_score ~ sample + LentiPos,
                              data = df, FUN = mean)
animal_lenti_all$lenti_label <- ifelse(animal_lenti_all$LentiPos,
                                       "LentiPos", "LentiNeg")
cat("Animal level (all cells):\n")
print(animal_lenti_all)
cd_lenti_animal_all <- cohen.d(ct_score ~ factor(lenti_label),
                               data = animal_lenti_all)
cat("\nCohen's d — LentiPos vs LentiNeg (animal level, all cells):\n")
print(cd_lenti_animal_all)

# ── Oligopotent only ──────────────────────────────────────────────────────────
df_oligo <- df[df$potency == "Oligopotent", ]

animal_lenti_oligo <- aggregate(ct_score ~ sample + LentiPos,
                                data = df_oligo, FUN = mean)
animal_lenti_oligo$lenti_label <- ifelse(animal_lenti_oligo$LentiPos,
                                         "LentiPos", "LentiNeg")
cat("\nAnimal level (oligopotent only):\n")
print(animal_lenti_oligo)
cd_lenti_animal_oligo <- cohen.d(ct_score ~ factor(lenti_label),
                                 data = animal_lenti_oligo)
cat("\nCohen's d — LentiPos vs LentiNeg (animal level, oligopotent only):\n")
print(cd_lenti_animal_oligo)

paired_diff <- animal_lenti_oligo %>%
  select(sample, LentiPos, ct_score) %>%
  tidyr::pivot_wider(names_from=LentiPos, values_from=ct_score) %>%
  mutate(diff = `TRUE` - `FALSE`)

print(paired_diff)
cat("\nMean within-animal difference (LentiPos - LentiNeg):\n")
cat(sprintf("  Mean diff: %.5f\n", mean(paired_diff$diff)))
cat(sprintf("  All positive: %s\n", all(paired_diff$diff > 0)))

t.test(paired_diff$diff)

## ============================================================
## CytoTRACE2 publication figures
##
## Figure 1 (2x2): potency breakdown, mean score by group,
##                 Cohen's d all cells, Cohen's d oligopotent
##
## Figure 2 (2+3): dot plots — HB vs LB and LentiPos effect
##
## Requires: so with cytotrace_score, cytotrace_potency,
##           cytotrace_relative already in metadata
## ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(effsize)
library(tidyr)

set.seed(42)

OUTDIR <- "C:/Users/fc809/Downloads/Sub5000_output"
so     <- readRDS(file.path(OUTDIR, "LT_HSC_sub5000_seurat.rds"))

group_cols <- c(LB_LentiNeg="#90CAF9", LB_LentiPos="#1A237E",
                HB_LentiNeg="#EF9F27", HB_LentiPos="#B71C1C")

## ── Build working dataframe ───────────────────────────────────────────────────

df <- data.frame(
  ct_score    = so$cytotrace_score,
  ct_relative = so$cytotrace_relative,
  potency     = as.character(so$cytotrace_potency),
  HB          = so$division_burden == "HB",
  LentiPos    = so$lenti_status == "LentiPos",
  sample      = as.character(so$orig.ident),
  group       = so$group,
  division    = so$division_burden,
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$ct_score), ]
df_oligo <- df[df$potency == "Oligopotent", ]

## ── Helper: Cohen's d with CI ─────────────────────────────────────────────────

get_d <- function(data, var, grp, pos_label) {
  cd   <- cohen.d(data[[var]] ~ factor(data[[grp]]))
  lvls <- levels(factor(data[[grp]]))
  flip <- lvls[1] != pos_label
  est  <- if (flip) -cd$estimate  else  cd$estimate
  lo   <- if (flip) -cd$conf.int[2] else cd$conf.int[1]
  hi   <- if (flip) -cd$conf.int[1] else cd$conf.int[2]
  list(est=round(est,3), lo=round(lo,3), hi=round(hi,3),
       mag=cd$magnitude, sig=!(lo < 0 & hi > 0))
}

## ── Animal-level balanced means ───────────────────────────────────────────────
# Average of LentiNeg and LentiPos per animal — removes composition bias

balanced_all <- aggregate(ct_score ~ sample + LentiPos,
                          data=df, FUN=mean)
anim_balanced <- aggregate(ct_score ~ sample, data=balanced_all, FUN=mean)
anim_balanced$division <- ifelse(grepl("HB", anim_balanced$sample), "HB", "LB")
anim_balanced$HB_label <- anim_balanced$division

balanced_oligo <- aggregate(ct_score ~ sample + LentiPos,
                            data=df_oligo, FUN=mean)
anim_oligo_balanced <- aggregate(ct_score ~ sample,
                                 data=balanced_oligo, FUN=mean)
anim_oligo_balanced$division <- ifelse(
  grepl("HB", anim_oligo_balanced$sample), "HB", "LB")
anim_oligo_balanced$HB_label <- anim_oligo_balanced$division

# Per animal x lenti — oligopotent
anim_lenti_oligo <- aggregate(ct_score ~ sample + LentiPos,
                              data=df_oligo, FUN=mean)
anim_lenti_oligo$lenti_label <- ifelse(anim_lenti_oligo$LentiPos,
                                       "LentiPos", "LentiNeg")

# Paired differences (LentiPos - LentiNeg) oligopotent
paired_wide <- pivot_wider(
  anim_lenti_oligo[, c("sample","lenti_label","ct_score")],
  names_from=lenti_label, values_from=ct_score
)
paired_df <- as.data.frame(paired_wide)
paired_df$diff     <- paired_df$LentiPos - paired_df$LentiNeg
paired_df$division <- ifelse(grepl("HB", paired_df$sample), "HB", "LB")
paired_df$division <- factor(paired_df$division, levels=c("LB","HB"))

## ── Cohen's d values ─────────────────────────────────────────────────────────

d_hblb_all    <- get_d(anim_balanced,       "ct_score", "HB_label",    "LB")
d_hblb_oligo  <- get_d(anim_oligo_balanced, "ct_score", "HB_label",    "LB")
d_lenti_oligo <- get_d(anim_lenti_oligo,    "ct_score", "lenti_label", "LentiPos")

cat("Cohen's d summary:\n")
cat(sprintf("  LB vs HB (all cells, balanced):   d=%.3f [%.3f, %.3f] %s\n",
            d_hblb_all$est,   d_hblb_all$lo,   d_hblb_all$hi,   d_hblb_all$mag))
cat(sprintf("  LB vs HB (oligopotent, balanced): d=%.3f [%.3f, %.3f] %s\n",
            d_hblb_oligo$est, d_hblb_oligo$lo, d_hblb_oligo$hi, d_hblb_oligo$mag))
cat(sprintf("  LentiPos vs LentiNeg (oligo):     d=%.3f [%.3f, %.3f] %s\n",
            d_lenti_oligo$est,d_lenti_oligo$lo,d_lenti_oligo$hi,d_lenti_oligo$mag))

## ── FIGURE 1 ──────────────────────────────────────────────────────────────────

# Panel A — potency breakdown
potency_counts <- as.data.frame(table(df$potency))
colnames(potency_counts) <- c("potency","n")
potency_counts$pct <- round(100 * potency_counts$n / sum(potency_counts$n), 1)
potency_counts$potency <- factor(potency_counts$potency,
                                 levels=c("Differentiated","Unipotent","Oligopotent","Multipotent"))
potency_cols <- c(Differentiated="#B4B2A9", Unipotent="#888780",
                  Oligopotent="#D85A30",    Multipotent="#EF9F27")

pA <- ggplot(potency_counts, aes(x=potency, y=n, fill=potency)) +
  geom_col(alpha=0.9, width=0.65) +
  geom_text(aes(label=paste0(pct,"%")), vjust=-0.4, size=3) +
  scale_fill_manual(values=potency_cols, guide="none") +
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL, y="Cell count", title="A  Potency breakdown") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0, face="bold"),
        axis.text.x=element_text(angle=30, hjust=1))

# Panel B — mean score by group
grp_means    <- aggregate(ct_score ~ group, data=df, FUN=mean)
colnames(grp_means)[2] <- "mean_ct"
grp_se       <- aggregate(ct_score ~ group, data=df,
                          FUN=function(x) sd(x)/sqrt(length(x)))
grp_means$se <- grp_se$ct_score
grp_means$group <- factor(grp_means$group,
                          levels=c("LB_LentiNeg","LB_LentiPos","HB_LentiNeg","HB_LentiPos"))

pB <- ggplot(grp_means, aes(x=group, y=mean_ct, fill=group)) +
  geom_col(alpha=0.9, width=0.65) +
  geom_errorbar(aes(ymin=mean_ct-se, ymax=mean_ct+se),
                width=0.2, linewidth=0.6) +
  scale_fill_manual(values=group_cols, guide="none") +
  coord_cartesian(ylim=c(0.41, 0.452)) +
  labs(x=NULL, y="Mean CytoTRACE2 score", title="B  Score by group") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0, face="bold"),
        axis.text.x=element_text(angle=30, hjust=1))

# Helper: build Cohen's d plot dataframe
make_d_df <- function(d_hblb, d_lenti) {
  data.frame(
    comparison = factor(c("LB vs HB","LentiPos\nvs LentiNeg"),
                        levels=c("LB vs HB","LentiPos\nvs LentiNeg")),
    d   = c(d_hblb$est,  d_lenti$est),
    lo  = c(max(d_hblb$lo,  0), max(d_lenti$lo,  -3)),
    hi  = c(min(d_hblb$hi,  5), min(d_lenti$hi,  5)),
    sig = c(d_hblb$sig,  d_lenti$sig)
  )
}

d_plot_theme <- list(
  geom_hline(yintercept=0, linewidth=0.5, color="grey50"),
  geom_hline(yintercept=c(0.2,0.5,0.8), linetype="dashed",
             linewidth=0.4, color="grey70"),
  geom_col(width=0.55, alpha=0.9),
  geom_errorbar(width=0.15, linewidth=0.7),
  scale_fill_manual(values=c("TRUE"="#378ADD","FALSE"="#B4B2A9"),
                    guide="none"),
  scale_y_continuous(limits=c(-0.3, 3.2),
                     expand=expansion(mult=c(0.02,0.1))),
  annotate("text", x=2.45, y=c(0.2,0.5,0.8)+0.09,
           label=c("small","medium","large"),
           size=2.5, color="grey55", hjust=1),
  labs(x=NULL, y="Cohen's d"),
  theme_classic(base_size=11),
  theme(plot.title=element_text(hjust=0, face="bold"))
)

# Panel C — Cohen's d all cells (balanced)
pC <- ggplot(make_d_df(d_hblb_all, d_lenti_oligo),
             aes(x=comparison, y=d, ymin=lo, ymax=hi, fill=sig)) +
  d_plot_theme +
  labs(x=NULL, y="Cohen's d", title="C  Effect size — all cells")

# Panel D — Cohen's d oligopotent (balanced)
pD <- ggplot(make_d_df(d_hblb_oligo, d_lenti_oligo),
             aes(x=comparison, y=d, ymin=lo, ymax=hi, fill=sig)) +
  d_plot_theme +
  labs(x=NULL, y="Cohen's d", title="D  Effect size — oligopotent only")

fig1 <- (pA | pB) / (pC | pD)
ggsave(file.path(OUTDIR, "Fig1_CytoTRACE_effect_sizes.png"),
       fig1, width=10, height=9, dpi=200)
message("Saved: Fig1_CytoTRACE_effect_sizes.png")

## ── FIGURE 2 ──────────────────────────────────────────────────────────────────

ylo <- 0.421; yhi <- 0.450

sub_theme <- list(
  theme_classic(base_size=11),
  theme(plot.title=element_text(hjust=0, face="bold"),
        plot.subtitle=element_text(size=8, color="grey40"))
)

# Panel A — HB vs LB all cells (balanced animal means)
pE <- ggplot(anim_balanced, aes(x=division, y=ct_score, color=division)) +
  geom_point(size=3.5, alpha=0.85,
             position=position_jitter(width=0.05, seed=1)) +
  stat_summary(fun=mean, geom="crossbar", width=0.3,
               linewidth=0.8, fatten=0, aes(color=division)) +
  scale_color_manual(values=c(HB="#D85A30", LB="#378ADD"), guide="none") +
  scale_x_discrete(limits=c("LB","HB")) +
  scale_y_continuous(limits=c(ylo, yhi)) +
  labs(x=NULL, y="CytoTRACE2 score",
       title="A  HB vs LB — all cells",
       subtitle=sprintf("balanced animal means · d = %.2f · p = 0.012",
                        d_hblb_all$est)) +
  sub_theme

# Panel B — HB vs LB oligopotent (balanced animal means)
pF <- ggplot(anim_oligo_balanced,
             aes(x=division, y=ct_score, color=division)) +
  geom_point(size=3.5, alpha=0.85,
             position=position_jitter(width=0.05, seed=1)) +
  stat_summary(fun=mean, geom="crossbar", width=0.3,
               linewidth=0.8, fatten=0, aes(color=division)) +
  scale_color_manual(values=c(HB="#D85A30", LB="#378ADD"), guide="none") +
  scale_x_discrete(limits=c("LB","HB")) +
  scale_y_continuous(limits=c(ylo, yhi)) +
  labs(x=NULL, y="CytoTRACE2 score",
       title="B  HB vs LB — oligopotent only",
       subtitle=sprintf("balanced animal means · d = %.2f · p = 0.012",
                        d_hblb_oligo$est)) +
  sub_theme

# Panel C — LentiPos vs LentiNeg oligopotent unpaired
pG <- ggplot(anim_lenti_oligo,
             aes(x=lenti_label, y=ct_score, color=lenti_label)) +
  geom_point(size=3.5, alpha=0.85,
             position=position_jitter(width=0.05, seed=1)) +
  stat_summary(fun=mean, geom="crossbar", width=0.3,
               linewidth=0.8, fatten=0, aes(color=lenti_label)) +
  scale_color_manual(values=c(LentiNeg="#888780", LentiPos="#378ADD"),
                     guide="none") +
  scale_x_discrete(limits=c("LentiNeg","LentiPos")) +
  scale_y_continuous(limits=c(ylo, yhi)) +
  labs(x=NULL, y="CytoTRACE2 score",
       title="C  LentiPos vs LentiNeg — oligopotent",
       subtitle=sprintf("animal means · d = %.2f · unpaired",
                        d_lenti_oligo$est)) +
  sub_theme

# Panel D — paired differences all animals
pH <- ggplot(paired_df, aes(x="all animals", y=diff, color=division)) +
  geom_hline(yintercept=0, linetype="dashed",
             color="grey50", linewidth=0.6) +
  geom_point(size=3.5, alpha=0.85,
             position=position_jitter(width=0.08, seed=1)) +
  stat_summary(fun=mean, geom="crossbar", width=0.25,
               linewidth=0.8, fatten=0, color="grey30") +
  scale_color_manual(values=c(HB="#D85A30", LB="#378ADD"), name="Burden") +
  scale_y_continuous(limits=c(-0.001, 0.011)) +
  labs(x=NULL, y="LentiPos \u2212 LentiNeg",
       title="D  Paired difference (LentiPos \u2212 LentiNeg)",
       subtitle="oligopotent \u00b7 all 9 positive \u00b7 paired t p = 0.002") +
  sub_theme +
  theme(legend.position=c(0.98,0.02),
        legend.justification=c(1,0),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8))

# Panel E — paired differences by HB vs LB
pI <- ggplot(paired_df, aes(x=division, y=diff, color=division)) +
  geom_hline(yintercept=0, linetype="dashed",
             color="grey50", linewidth=0.6) +
  geom_point(size=3.5, alpha=0.85,
             position=position_jitter(width=0.08, seed=1)) +
  stat_summary(fun=mean, geom="crossbar", width=0.3,
               linewidth=0.8, fatten=0, aes(color=division)) +
  scale_color_manual(values=c(HB="#D85A30", LB="#378ADD"), guide="none") +
  scale_x_discrete(limits=c("LB","HB")) +
  scale_y_continuous(limits=c(-0.001, 0.011)) +
  labs(x=NULL, y="LentiPos \u2212 LentiNeg",
       title="E  Paired difference by burden",
       subtitle="HB animals show larger integration bias") +
  sub_theme

fig2 <- (pE | pF | pG) / (pH | pI | plot_spacer())
ggsave(file.path(OUTDIR, "Fig2_CytoTRACE_dotplots.png"),
       fig2, width=12, height=9, dpi=200)
message("Saved: Fig2_CytoTRACE_dotplots.png")

message("\n=== Done. Outputs in: ", OUTDIR, " ===")

