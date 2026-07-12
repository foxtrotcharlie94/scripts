# ============================================================
# BMDM data — full diagnostic + concordance analysis (one script)
#
# Works entirely from the files already in the BMDM folder:
#   <comp>/DEG_<comp>_bulk.txt                         (edgeR results + per-sample CPM)
#   <comp>/edgeR_<comp>_TMM_normalized_CPM_wOffset.csv (TMM-normalized CPM)
# No raw counts required.
#
# Four comparisons (folder -> meaning; sign of logFC):
#   HB.neg_vs_LB.neg  : HB vs LB within lenti-neg   (+ = up in HB)
#   HB.pos_vs_LB.pos  : HB vs LB within lenti-pos   (+ = up in HB)
#   HB.pos_vs_HB.neg  : pos vs neg within HB        (+ = up in lenti-pos)
#   LB.pos_vs_LB.neg  : pos vs neg within LB        (+ = up in lenti-pos)
#
# Sections:
#   1. DGE asymmetry summary (up/down at FDR thresholds)
#   2. logFC distribution + tail asymmetry (histograms + ranked waterfall)
#   3. Library-composition check (CPM-based) + abundant-gene leverage
#   4. Intensity dependence: logFC vs logCPM
#   5. Cross-axis concordance: burden (HB vs LB) vs lenti (pos vs neg)
#        5a. all genes        (Pearson, Spearman, R^2, regression)
#        5b. |log2FC|>0.5 both (same + quadrant %)
#   6. CH / myeloid-inflammation gene heatmap across the 4 groups
#        (expression level + paired/unpaired t-tests, BH-corrected)
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  have_repel <- requireNamespace("ggrepel", quietly = TRUE)
})
theme_set(theme_bw(base_size = 12))

# ---- settings ----
bmdm_dir <- "C:/Users/fc809/Downloads/BMDM"
out_dir  <- file.path(bmdm_dir, "axis_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COMPS <- c("HB.neg_vs_LB.neg","HB.pos_vs_LB.pos","HB.pos_vs_HB.neg","LB.pos_vs_LB.neg")
POSMEANS <- c(HB.neg_vs_LB.neg="+ = up in HB",  HB.pos_vs_LB.pos="+ = up in HB",
              HB.pos_vs_HB.neg="+ = up in pos", LB.pos_vs_LB.neg="+ = up in pos")
F_NEG <- "HB.neg_vs_LB.neg"; F_POS <- "HB.pos_vs_LB.pos"   # group-mean sources for heatmap/composition

# ---- readers ----
read_deg <- function(f) {
  d <- fread(file.path(bmdm_dir, f, sprintf("DEG_%s_bulk.txt", f)))
  setnames(d, 1, "ensembl")
  d[, .(ensembl, gene = gene_name, logFC, logCPM, F, PValue, FDR)]
}
read_cpm <- function(f) {
  d <- fread(file.path(bmdm_dir, f, sprintf("edgeR_%s_TMM_normalized_CPM_wOffset.csv", f)))
  setnames(d, 1, "ensembl")
  m <- as.matrix(d[, -1]); rownames(m) <- d$ensembl; m
}
deg <- setNames(lapply(COMPS, read_deg), COMPS)
cpm <- setNames(lapply(COMPS, read_cpm), COMPS)
grepcols <- function(m, pat) grep(pat, colnames(m), value = TRUE)

# ════════════════════════════════════════════════════════════════════════════
# 1. DGE ASYMMETRY SUMMARY
# ════════════════════════════════════════════════════════════════════════════
sec1 <- rbindlist(lapply(COMPS, function(f) {
  d <- deg[[f]]
  s05 <- d[FDR < 0.05]; s25 <- d[FDR < 0.25]
  data.table(comparison = f, sign = POSMEANS[[f]],
             n_genes = nrow(d), min_FDR = signif(min(d$FDR, na.rm = TRUE), 3),
             up_FDR05 = sum(s05$logFC > 0), dn_FDR05 = sum(s05$logFC < 0),
             up_FDR25 = sum(s25$logFC > 0), dn_FDR25 = sum(s25$logFC < 0),
             med_logFC = round(median(d$logFC), 3), mean_logFC = round(mean(d$logFC), 3))
}))
fwrite(sec1, file.path(out_dir, "1_DGE_asymmetry_summary.csv"))
cat("== 1. DGE asymmetry ==\n"); print(sec1)

# bar plot of up/down at FDR<0.25
bar <- melt(sec1[, .(comparison, up = up_FDR25, down = dn_FDR25)],
            id.vars = "comparison", variable.name = "dir", value.name = "n")
bar[, comparison := factor(comparison, levels = COMPS)]
p1 <- ggplot(bar, aes(comparison, ifelse(dir == "down", -n, n), fill = dir)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_fill_manual(values = c(up = "#C0392B", down = "#2C6FBF"), name = NULL) +
  labs(title = "DGE asymmetry — significant genes (FDR < 0.25)",
       x = NULL, y = "n genes (up above / down below 0)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out_dir, "1_DGE_asymmetry.png"), p1, width = 8, height = 5, dpi = 150)

# ════════════════════════════════════════════════════════════════════════════
# 2. logFC DISTRIBUTION + TAIL ASYMMETRY
# ════════════════════════════════════════════════════════════════════════════
allfc <- rbindlist(lapply(COMPS, function(f)
  data.table(comparison = f, gene = deg[[f]]$gene, logFC = deg[[f]]$logFC)))
allfc[, comparison := factor(comparison, levels = COMPS)]

# top-N |logFC| up/down split
sec2 <- rbindlist(lapply(COMPS, function(f) {
  d <- deg[[f]][order(-abs(logFC))]
  tn <- function(n) c(up = sum(head(d$logFC, n) > 0), dn = sum(head(d$logFC, n) < 0))
  t50 <- tn(50); t200 <- tn(200)
  data.table(comparison = f, top50_up = t50["up"], top50_dn = t50["dn"],
             top200_up = t200["up"], top200_dn = t200["dn"])
}))
fwrite(sec2, file.path(out_dir, "2_logFC_tail_split.csv"))
cat("\n== 2. tail split ==\n"); print(sec2)

# histograms
p2a <- ggplot(allfc, aes(logFC)) +
  geom_histogram(bins = 80, fill = "grey60") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  facet_wrap(~ comparison, scales = "free_y", nrow = 1) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(title = "logFC distribution per comparison", x = "logFC", y = "# genes")
# ranked waterfall
wf <- allfc[, .(logFC = sort(logFC, decreasing = TRUE)), by = comparison]
wf[, rank := seq_len(.N), by = comparison]
p2b <- ggplot(wf, aes(rank, logFC, fill = logFC > 0)) +
  geom_col(width = 1) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2C6FBF"), guide = "none") +
  facet_wrap(~ comparison, scales = "free", nrow = 1) +
  labs(title = "Ranked logFC (waterfall) — tail asymmetry", x = "gene rank", y = "logFC (sorted)")
ggsave(file.path(out_dir, "2_logFC_distribution.png"), p2a / p2b,
       width = 16, height = 8, dpi = 150)

# ════════════════════════════════════════════════════════════════════════════
# 3. LIBRARY-COMPOSITION CHECK  (CPM/1e6 = fraction of library)
# ════════════════════════════════════════════════════════════════════════════
# assemble one column per sample: neg samples from F_NEG, pos samples from F_POS
samp_mat <- cbind(cpm[[F_NEG]][, grepcols(cpm[[F_NEG]], "neg$"), drop = FALSE],
                  cpm[[F_POS]][, grepcols(cpm[[F_POS]], "pos$"), drop = FALSE])
group_of <- function(s) paste0(ifelse(grepl("^HB", s), "HB", "LB"),
                               ".lenti", ifelse(grepl("pos$", s), "pos", "neg"))
GROUPS <- c("HB.lentineg","HB.lentipos","LB.lentineg","LB.lentipos")

comp_tbl <- rbindlist(lapply(colnames(samp_mat), function(s) {
  v <- sort(samp_mat[, s], decreasing = TRUE) / 1e6
  data.table(sample = s, group = group_of(s),
             top1 = v[1], top5 = sum(v[1:5]), top10 = sum(v[1:10]),
             top20 = sum(v[1:20]), top50 = sum(v[1:50]),
             top1_gene = deg[[F_NEG]]$gene[match(names(v)[1], deg[[F_NEG]]$ensembl)])
}))
fwrite(comp_tbl, file.path(out_dir, "3_composition_per_sample.csv"))
cat("\n== 3. composition (group means, % of library) ==\n")
print(comp_tbl[, lapply(.SD, function(x) round(mean(x)*100,1)),
               by = group, .SDcols = c("top1","top5","top10","top20","top50")])

# dominance plot per sample
cm <- melt(comp_tbl, id.vars = c("sample","group"),
           measure.vars = c("top1","top5","top10","top20","top50"),
           variable.name = "level", value.name = "frac")
cm[, sample := factor(sample, levels = comp_tbl[order(match(group, GROUPS)), sample])]
p3a <- ggplot(cm, aes(sample, 100*frac, color = level, group = level)) +
  geom_line() + geom_point(size = 1.5) +
  labs(title = "Library dominance by top-N genes (per sample)",
       x = NULL, y = "% of library (TMM-normalized CPM)", color = "cumulative") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
# cumulative share (group mean)
cum <- rbindlist(lapply(GROUPS, function(g) {
  cols <- comp_tbl[group == g, sample]
  v <- sort(rowMeans(samp_mat[, cols, drop = FALSE]), decreasing = TRUE) / 1e6
  data.table(group = g, n = seq_along(v), cum = cumsum(v)*100)
}))
p3b <- ggplot(cum[n <= 2000], aes(n, cum, color = group)) +
  geom_line(linewidth = 0.9) + scale_x_log10() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
  labs(title = "Cumulative library share (group mean)",
       x = "number of top genes (log)", y = "cumulative % of library")
ggsave(file.path(out_dir, "3_composition_check.png"), p3a | p3b, width = 14, height = 5.5, dpi = 150)

# 3b. abundant-gene leverage: are the most abundant genes up in the numerator,
#     and do they take a larger library share there?
sec3b <- rbindlist(lapply(COMPS, function(f) {
  m <- cpm[[f]]; d <- deg[[f]]
  parts <- strsplit(f, "_vs_")[[1]]                  # e.g. c("HB.pos","HB.neg")
  half  <- function(h) { cond <- substr(h, 1, 2); geno <- sub("^(HB|LB)\\.", "", h)
                         grepcols(m, paste0("^", cond, "[0-9]+", geno, "$")) }
  A <- half(parts[1]); B <- half(parts[2])           # A = numerator group
  mA <- rowMeans(m[, A, drop = FALSE]); mB <- rowMeans(m[, B, drop = FALSE])
  ord <- names(sort((mA + mB) / 2, decreasing = TRUE))      # genes by overall abundance
  lfc <- d$logFC[match(ord, d$ensembl)]
  shareA <- function(n) sum(sort(mA, decreasing = TRUE)[1:n]) / 1e6
  shareB <- function(n) sum(sort(mB, decreasing = TRUE)[1:n]) / 1e6
  data.table(comparison = f,
             med_logFC_top20abund  = round(median(lfc[1:20],  na.rm = TRUE), 3),
             med_logFC_top100abund = round(median(lfc[1:100], na.rm = TRUE), 3),
             top50_share_A = round(shareA(50)*100, 1), top50_share_B = round(shareB(50)*100, 1),
             top50_share_diff_pp = round((shareA(50) - shareB(50))*100, 1))
}))
fwrite(sec3b, file.path(out_dir, "3b_abundant_gene_leverage.csv"))
cat("\n== 3b. abundant-gene leverage ==\n"); print(sec3b)

# ════════════════════════════════════════════════════════════════════════════
# 4. INTENSITY DEPENDENCE: logFC vs logCPM
# ════════════════════════════════════════════════════════════════════════════
sec4 <- rbindlist(lapply(COMPS, function(f) {
  d <- deg[[f]][is.finite(logFC) & is.finite(logCPM)]
  ter <- cut(d$logCPM, quantile(d$logCPM, c(0,1/3,2/3,1)), labels = c("low","mid","high"),
             include.lowest = TRUE)
  mt <- tapply(d$logFC, ter, mean)
  data.table(comparison = f,
             pearson = round(cor(d$logCPM, d$logFC), 3),
             spearman = round(cor(d$logCPM, d$logFC, method = "spearman"), 3),
             meanlogFC_low = round(mt["low"],3), meanlogFC_mid = round(mt["mid"],3),
             meanlogFC_high = round(mt["high"],3))
}))
fwrite(sec4, file.path(out_dir, "4_intensity_dependence.csv"))
cat("\n== 4. intensity dependence (logFC vs logCPM) ==\n"); print(sec4)

# intensity plot (logFC vs logCPM) faceted
int_dt <- rbindlist(lapply(COMPS, function(f)
  data.table(comparison = f, logCPM = deg[[f]]$logCPM, logFC = deg[[f]]$logFC)))
int_dt[, comparison := factor(comparison, levels = COMPS)]
p4 <- ggplot(int_dt, aes(logCPM, logFC)) +
  geom_point(size = 0.3, alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, color = "#C0392B", linewidth = 0.8) +
  facet_wrap(~ comparison, nrow = 1) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(title = "Intensity dependence: logFC vs expression (logCPM)",
       subtitle = "downward trend at low expression = intensity-dependent bias")
ggsave(file.path(out_dir, "4_intensity_dependence.png"), p4, width = 16, height = 4.5, dpi = 150)

# ════════════════════════════════════════════════════════════════════════════
# 5. CROSS-AXIS CONCORDANCE
#    burden = avg logFC of (HB.neg_vs_LB.neg, HB.pos_vs_LB.pos)   + = up in HB
#    lenti  = avg logFC of (HB.pos_vs_HB.neg, LB.pos_vs_LB.neg)   + = up in lenti-pos
# ════════════════════════════════════════════════════════════════════════════
m <- Reduce(function(a,b) merge(a,b,by="ensembl"),
            list(deg[["HB.neg_vs_LB.neg"]][,.(ensembl,gene,b1=logFC)],
                 deg[["HB.pos_vs_LB.pos"]][,.(ensembl,b2=logFC)],
                 deg[["HB.pos_vs_HB.neg"]][,.(ensembl,l1=logFC)],
                 deg[["LB.pos_vs_LB.neg"]][,.(ensembl,l2=logFC)]))
m[, burden := (b1+b2)/2][, lenti := (l1+l2)/2]
m <- m[is.finite(burden) & is.finite(lenti)]
fwrite(m[,.(ensembl,gene,burden,lenti)], file.path(out_dir, "5_axis_logFC.csv"))

XL <- "Burden effect: HB vs LB  (log2FC, + = up in HB)"
YL <- "Lenti effect: pos vs neg  (log2FC, + = up in lentipos)"
ann <- function(r, rho, n, extra="") sprintf("n = %d\nPearson r = %.2f\nSpearman = %.2f\nR^2 = %.2f%s",
                                             n, r, rho, r^2, extra)

# 5a. all genes
r <- cor(m$burden, m$lenti); rho <- cor(m$burden, m$lenti, method = "spearman")
fit <- coef(lm(lenti ~ burden, m))
lim <- ceiling(quantile(abs(c(m$burden, m$lenti)), 0.995)*2)/2
p5a <- ggplot(m, aes(burden, lenti)) +
  geom_point(size = 0.4, alpha = 0.10, color = "#34618E") +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey65") +
  geom_abline(slope = fit[2], intercept = fit[1], color = "#C0392B", linewidth = 1) +
  annotate("label", x = -lim, y = lim, hjust = 0, vjust = 1, label = ann(r, rho, nrow(m)), size = 4) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Axis concordance — all genes", x = XL, y = YL)
ggsave(file.path(out_dir, "5a_concordance_allgenes.png"), p5a, width = 7.5, height = 7.5, dpi = 160)

# 5b. movers |log2FC|>0.5 on BOTH axes, with quadrant %
TH <- 0.5
mv <- m[abs(burden) > TH & abs(lenti) > TH]
rM <- cor(mv$burden, mv$lenti); rhoM <- cor(mv$burden, mv$lenti, method = "spearman")
fitM <- coef(lm(lenti ~ burden, mv)); N <- nrow(mv)
mv[, quad := fifelse(burden>0 & lenti>0, "up_up",
              fifelse(burden<0 & lenti<0, "dn_dn", "disc"))]
qp <- function(q) round(100*sum(mv$quad==q)/N)
conc <- qp("up_up") + qp("dn_dn")
mv[, col := fifelse(quad=="up_up","#C0392B", fifelse(quad=="dn_dn","#2C6FBF","#9AA0A6"))]
lim2 <- max(ceiling(quantile(abs(c(mv$burden,mv$lenti)),0.99)*2)/2, 1.5)
lab <- mv[quad %in% c("up_up","dn_dn")][order(-(abs(burden)+abs(lenti)))][1:min(18,.N)]
p5b <- ggplot(mv, aes(burden, lenti)) +
  annotate("rect", xmin=0, xmax=lim2, ymin=0, ymax=lim2, fill="#C0392B", alpha=0.06) +
  annotate("rect", xmin=-lim2, xmax=0, ymin=-lim2, ymax=0, fill="#2C6FBF", alpha=0.06) +
  geom_hline(yintercept=0,color="grey60") + geom_vline(xintercept=0,color="grey60") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey65") +
  geom_point(aes(color = col), size = 2) + scale_color_identity() +
  geom_abline(slope=fitM[2], intercept=fitM[1], color="black", linewidth=1) +
  annotate("label", x=-lim2, y=lim2, hjust=0, vjust=1, size = 4,
           label = ann(rM, rhoM, N,
                       sprintf("\nsame direction = %d%%\ndown/down %d%% (n=%d)\nup/up %d%% (n=%d)",
                               conc, qp("dn_dn"), sum(mv$quad=="dn_dn"),
                               qp("up_up"), sum(mv$quad=="up_up")))) +
  coord_equal(xlim=c(-lim2,lim2), ylim=c(-lim2,lim2)) +
  labs(title = sprintf("Axis concordance — responding genes (|log2FC|>%.1f on both)", TH),
       x = XL, y = YL)
if (have_repel)
  p5b <- p5b + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, fontface = "italic", max.overlaps = 30)
ggsave(file.path(out_dir, "5b_concordance_movers.png"), p5b, width = 7.5, height = 7.5, dpi = 160)
cat(sprintf("\n== 5. concordance ==\n  all: r=%.2f rho=%.2f R2=%.2f n=%d\n  movers: r=%.2f rho=%.2f R2=%.2f n=%d | up/up=%d%% dn/dn=%d%% concordant=%d%%\n",
            r, rho, r^2, nrow(m), rM, rhoM, rM^2, N, qp("up_up"), qp("dn_dn"), conc))

# ════════════════════════════════════════════════════════════════════════════
# 6. CH / MYELOID-INFLAMMATION GENE HEATMAP (4 groups) + paired/unpaired t-tests
# ════════════════════════════════════════════════════════════════════════════
GENES <- c("Il1b","Il1a","Il6","Tnf","Il18","Nlrp3","Casp1","Pycard","Aim2",
           "Ccl2","Ccl5","Cxcl1","Cxcl2","Cxcl10","Tlr2","Tlr4","Nfkb1","Tnfaip3","Ptgs2","Nos2")
sym2ens <- deg[[F_NEG]][, setNames(ensembl, gene)]
for (f in COMPS) { mp <- deg[[f]][, setNames(ensembl, gene)]
  sym2ens <- c(sym2ens, mp[setdiff(names(mp), names(sym2ens))]) }
ens <- sym2ens[GENES]

# group mean log2(CPM+1): neg groups from F_NEG, pos groups from F_POS
gmean <- function(cpmfile, pat) {
  m <- cpm[[cpmfile]]; cols <- grepcols(m, pat)
  v <- setNames(rep(NA_real_, length(GENES)), GENES)
  ok <- !is.na(ens) & ens %in% rownames(m)
  v[ok] <- rowMeans(log2(m[ens[ok], cols, drop = FALSE] + 1))
  v
}
gmean_lin <- function(cpmfile, pat) {                 # mean LINEAR CPM (for cell labels)
  m <- cpm[[cpmfile]]; cols <- grepcols(m, pat)
  v <- setNames(rep(NA_real_, length(GENES)), GENES)
  ok <- !is.na(ens) & ens %in% rownames(m)
  v[ok] <- rowMeans(m[ens[ok], cols, drop = FALSE])
  v
}
L <- cbind(HB.lentineg = gmean(F_NEG, "^HB[0-9]+neg$"),
           HB.lentipos = gmean(F_POS, "^HB[0-9]+pos$"),
           LB.lentineg = gmean(F_NEG, "^LB[0-9]+neg$"),
           LB.lentipos = gmean(F_POS, "^LB[0-9]+pos$"))
rownames(L) <- GENES
Clin <- cbind(HB.lentineg = gmean_lin(F_NEG, "^HB[0-9]+neg$"),
              HB.lentipos = gmean_lin(F_POS, "^HB[0-9]+pos$"),
              LB.lentineg = gmean_lin(F_NEG, "^LB[0-9]+neg$"),
              LB.lentipos = gmean_lin(F_POS, "^LB[0-9]+pos$"))
rownames(Clin) <- GENES
Z <- t(scale(t(L)))                                   # row z-score for color

# 4 contrasts (each from its co-normalized file). + = up in HB / up in pos
CT <- list(
  list(lab="HB vs LB | neg",  file="HB.neg_vs_LB.neg", A="^HB[0-9]+neg$", B="^LB[0-9]+neg$", paired=FALSE),
  list(lab="HB vs LB | pos",  file="HB.pos_vs_LB.pos", A="^HB[0-9]+pos$", B="^LB[0-9]+pos$", paired=FALSE),
  list(lab="pos vs neg | HB", file="HB.pos_vs_HB.neg", A="^HB[0-9]+pos$", B="^HB[0-9]+neg$", paired=TRUE),
  list(lab="pos vs neg | LB", file="LB.pos_vs_LB.neg", A="^LB[0-9]+pos$", B="^LB[0-9]+neg$", paired=TRUE))
donor <- function(s) sub("(pos|neg)$","",s)
eff <- matrix(NA, length(GENES), length(CT), dimnames=list(GENES, sapply(CT,`[[`,"lab")))
pv  <- eff
for (j in seq_along(CT)) {
  ct <- CT[[j]]; m <- cpm[[ct$file]]
  A <- grepcols(m, ct$A); B <- grepcols(m, ct$B)
  if (ct$paired) { da <- setNames(A, donor(A)); db <- setNames(B, donor(B))
                   d <- intersect(names(da), names(db)); A <- da[d]; B <- db[d] }
  for (i in seq_along(GENES)) {
    if (is.na(ens[i]) || !(ens[i] %in% rownames(m))) next   # leave NA if gene absent
    va <- log2(m[ens[i], A] + 1); vb <- log2(m[ens[i], B] + 1)
    eff[i,j] <- mean(va) - mean(vb)
    pv[i,j]  <- tryCatch(t.test(va, vb, paired = ct$paired)$p.value, error=function(e) NA)
  }
}
qv <- apply(pv, 2, p.adjust, method = "BH")           # BH across 20 genes per contrast
stars <- function(q) ifelse(is.na(q),"",ifelse(q<.001,"***",ifelse(q<.01,"**",ifelse(q<.05,"*",ifelse(q<.10,"·","")))))

ord <- GENES                                           # keep curated order
exp_dt <- as.data.table(as.table(Z)); setnames(exp_dt, c("gene","group","z"))
exp_dt[, gene := factor(gene, levels = rev(ord))]
exp_dt[, group := factor(group, levels = colnames(L))]
exp_dt[, cpm := Clin[cbind(as.character(gene), as.character(group))]]      # mean linear CPM
exp_dt[, cpm_lab := fifelse(is.na(cpm), "", fifelse(cpm >= 10, sprintf("%.0f", cpm), sprintf("%.1f", cpm)))]
exp_dt[, txt := fifelse(!is.na(z) & abs(z) > 1.2, "white", "black")]
pE <- ggplot(exp_dt, aes(group, gene, fill = z)) +
  geom_tile(color = "white") +
  geom_text(aes(label = cpm_lab, color = txt), size = 2.7) +
  scale_color_identity() +
  scale_fill_gradient2(low="#2C6FBF", mid="white", high="#C0392B", midpoint=0, name="z-score") +
  labs(title="Expression level — numbers = mean CPM\n(fill = row z-score of mean log2 CPM)", x=NULL, y=NULL) +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle=30, hjust=1),
                                        axis.text.y = element_text(face="italic"))

eff_dt <- as.data.table(as.table(eff)); setnames(eff_dt, c("gene","contrast","effect"))
eff_dt[, star := stars(as.vector(qv))]
eff_dt[, gene := factor(gene, levels = rev(ord))]
eff_dt[, contrast := factor(contrast, levels = colnames(eff))]
eff_dt[, lab := fifelse(is.na(effect), "", paste0(sprintf("%.1f", effect), star))]    # log2 diff + stars
eff_dt[, txt := fifelse(!is.na(effect) & abs(effect) > 1.0, "white", "black")]
pS <- ggplot(eff_dt, aes(contrast, gene, fill = effect)) +
  geom_tile(color = "white") +
  geom_text(aes(label = lab, color = txt), size = 2.7) +
  scale_color_identity() +
  geom_vline(xintercept = 2.5, color = "grey25", linewidth = 0.7) +
  