# ============================================================================
# Clock CpG correlation analysis
#
# Question: does forcing HSCs to proliferate (low-dose transplant +
# reconstitution) recapitulate, at the level of individual CpGs, the same
# methylation changes seen with natural chronological aging?
#
# This reproduces the reanalysis behind Fig 2E of Gorelov et al. 2024
# (Stem Cell Reports), using the raw GEO dataset it draws on (GSE44117,
# Beerman et al. 2013) and the clock CpG coefficient tables published in
# Meer et al. 2018 (eLife), for three independently-trained mouse DNAm
# clocks: Blood (Petkovich 2017), WLMT (Meer 2018), and YOMT (Stubbs 2017).
#
# See REPORT.md in this folder for full methods, data provenance, and
# results. Re-running this script regenerates every file in this folder.
# ============================================================================

## ---- 0. Paths & packages -------------------------------------------------

# Edit this if your folder is somewhere else.
base_dir    <- "C:/Users/fc809/Downloads/division_burden_methylation_analysis"
out_dir     <- file.path(base_dir, "clock_correlation_analysis")
extract_dir <- file.path(base_dir, "extracted")

tar_path    <- file.path(base_dir, "GSE44117_RAW.tar")
xlsx_path   <- file.path(base_dir, "elife-40675-supp3-v2.xlsx")
chain_path  <- file.path(base_dir, "mm10ToMm9.over.chain.gz")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_cran <- c("readxl", "data.table", "ggplot2", "BiocManager")
required_bioc <- c("rtracklayer", "GenomicRanges")

missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

library(readxl)
library(data.table)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)

## ---- 1. Extract raw RRBS data (GSE44117), if not already done -----------

if (!dir.exists(extract_dir)) dir.create(extract_dir, recursive = TRUE)
if (length(list.files(extract_dir, pattern = "\\.bed\\.gz$")) < 25) {
  message("Extracting ", tar_path, " ...")
  untar(tar_path, exdir = extract_dir)
}

# The 7 samples relevant to Fig 2E: young/old baseline HSCs, and HSCs from
# young donors transplanted at low dose (10 cells) and reconstituted for
# 20 weeks (forced proliferation).
samples <- c(
  young_1   = "GSM1079935_RRBS_cpgMethylation_Mouse_blood_HSC_young_1.RRBS.bed.gz",
  young_2   = "GSM1079939_RRBS_cpgMethylation_Mouse_blood_HSC_young_2.RRBS.bed.gz",
  old_3     = "GSM1079926_RRBS_cpgMethylation_Mouse_blood_HSC_old_3.RRBS.bed.gz",
  old_4     = "GSM1079927_RRBS_cpgMethylation_Mouse_blood_HSC_old_4.RRBS.bed.gz",
  reconst_1 = "GSM1079936_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_1.RRBS.bed.gz",
  reconst_2 = "GSM1079937_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_2.RRBS.bed.gz",
  reconst_3 = "GSM1079938_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_3.RRBS.bed.gz"
)
young_samples   <- c("young_1", "young_2")
old_samples     <- c("old_3", "old_4")
reconst_samples <- c("reconst_1", "reconst_2", "reconst_3")

## ---- 2. Load clock CpG coefficient tables (mm10) -------------------------

blood_tab <- as.data.frame(read_excel(xlsx_path, sheet = "Blood"))
wlmt_tab  <- as.data.frame(read_excel(xlsx_path, sheet = "Whole lifespan multi-tissue"))
yomt_tab  <- as.data.frame(read_excel(xlsx_path, sheet = "Young age multi-tissue"))
# (The "Liver" sheet is also in this file but is skipped: not relevant to
#  blood/HSC samples.)

tidy_cpg <- function(tab) {
  d <- data.frame(chr = tab$Chromosome, pos = as.integer(tab$Position),
                   weight = as.numeric(tab$Weight))
  d[!is.na(d$chr) & !is.na(d$pos), ]
}
blood_cpg <- tidy_cpg(blood_tab)
wlmt_cpg  <- tidy_cpg(wlmt_tab)
yomt_cpg  <- tidy_cpg(yomt_tab)

## ---- 3. LiftOver clock coordinates mm10 -> mm9 ---------------------------
# GSE44117 was mapped to mm9 (Maq aligner, standard for that 2012-13 study);
# the clock tables are on mm10. This was confirmed empirically: before
# liftover, ~0% of clock CpGs matched RRBS-covered positions; after, 75-93%
# did.

chain_path_plain <- sub("\\.gz$", "", chain_path)
if (!file.exists(chain_path_plain)) {
  # import.chain() does not auto-decompress .gz files (it reads the raw
  # gzip bytes as text, causing "expected 11 elements in header, got 1").
  # Decompress once to a plain-text .chain file via gzfile(), which does
  # handle gzip properly, and reuse it on subsequent runs.
  gz_con <- gzfile(chain_path, "rt")
  chain_lines <- readLines(gz_con)
  close(gz_con)
  writeLines(chain_lines, chain_path_plain)
}
chain <- import.chain(chain_path_plain)

liftover_cpgs <- function(cpg_df) {
  # Excel "Position" is 1-based; liftOver here on the 0-based coordinate,
  # matching how this was validated against the RRBS BED coordinates.
  gr <- GRanges(seqnames = cpg_df$chr, ranges = IRanges(start = cpg_df$pos - 1, width = 1))
  lifted <- liftOver(gr, chain)
  keep <- lengths(lifted) == 1          # drop CpGs with no/ambiguous mapping
  lifted_kept <- unlist(lifted[keep])   # order-preserving since each has exactly 1 range
  out <- cpg_df[keep, ]
  out$pos_mm9 <- start(lifted_kept)
  out
}

blood_mm9 <- liftover_cpgs(blood_cpg)
wlmt_mm9  <- liftover_cpgs(wlmt_cpg)
yomt_mm9  <- liftover_cpgs(yomt_cpg)

message(sprintf("Lifted to mm9: Blood %d/%d, WLMT %d/%d, YOMT %d/%d",
                 nrow(blood_mm9), nrow(blood_cpg),
                 nrow(wlmt_mm9), nrow(wlmt_cpg),
                 nrow(yomt_mm9), nrow(yomt_cpg)))

## ---- 4. Parse RRBS BED files & match clock CpGs to methylation values ----
# RRBS reports one row per covered CpG as chr/start/end/'methylated/total'/
# score/strand. Because which strand's read covers a given CpG can shift the
# reported coordinate by 1 bp, we try offset -1 then fall back to offset 0.

parse_sample <- function(path, cpg_mm9) {
  dt <- fread(path, header = FALSE,
              col.names = c("chr", "start", "end", "ratio", "score", "strand"))
  dt[, c("meth", "total") := tstrsplit(gsub("'", "", ratio), "/", fixed = TRUE)]
  dt[, meth := as.numeric(meth)]
  dt[, total := as.numeric(total)]
  dt[, pct := 100 * meth / total]

  key1 <- data.table(chr = cpg_mm9$chr, start = cpg_mm9$pos_mm9 - 1)
  key0 <- data.table(chr = cpg_mm9$chr, start = cpg_mm9$pos_mm9)

  m1 <- dt[key1, on = c("chr", "start")]$pct
  m0 <- dt[key0, on = c("chr", "start")]$pct
  ifelse(!is.na(m1), m1, m0)
}

compute_deltas <- function(cpg_mm9, clock_name) {
  vals <- sapply(names(samples), function(s) parse_sample(file.path(extract_dir, samples[[s]]), cpg_mm9))
  colnames(vals) <- names(samples)

  keep <- complete.cases(vals)
  vals <- vals[keep, , drop = FALSE]
  cpg_kept <- cpg_mm9[keep, ]

  young_mean   <- rowMeans(vals[, young_samples, drop = FALSE])
  old_mean     <- rowMeans(vals[, old_samples, drop = FALSE])
  reconst_mean <- rowMeans(vals[, reconst_samples, drop = FALSE])

  data.frame(
    clock = clock_name, chr = cpg_kept$chr, pos = cpg_kept$pos_mm9, weight = cpg_kept$weight,
    delta_age_old_minus_young = old_mean - young_mean,
    delta_reconst_minus_young = reconst_mean - young_mean
  )
}

message("Parsing RRBS samples and computing per-CpG deltas (may take a minute)...")
df_all <- rbind(
  compute_deltas(blood_mm9, "Blood"),
  compute_deltas(wlmt_mm9,  "WLMT"),
  compute_deltas(yomt_mm9,  "YOMT")
)
fwrite(df_all, file.path(out_dir, "cpg_deltas.csv"))
message(sprintf("CpGs with full coverage across all 7 samples: Blood=%d, WLMT=%d, YOMT=%d",
                 sum(df_all$clock == "Blood"), sum(df_all$clock == "WLMT"), sum(df_all$clock == "YOMT")))

## ---- 5. Correlation & regression, per clock -------------------------------

fit_stats <- do.call(rbind, lapply(split(df_all, df_all$clock), function(sub) {
  fit <- lm(delta_reconst_minus_young ~ delta_age_old_minus_young, data = sub)
  ct  <- cor.test(sub$delta_age_old_minus_young, sub$delta_reconst_minus_young)
  data.frame(clock = sub$clock[1], n = nrow(sub),
             slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]),
             r = unname(ct$estimate), p = ct$p.value)
}))
print(fit_stats)
fwrite(fit_stats, file.path(out_dir, "fit_stats.csv"))

## ---- 6. Plots --------------------------------------------------------------

clock_colors <- c(Blood = "#B3B3B3", WLMT = "#1B3A5C", YOMT = "#CC9900")

fit_line <- function(sub, xrange) {
  fit <- lm(delta_reconst_minus_young ~ delta_age_old_minus_young, data = sub)
  xs <- seq(xrange[1], xrange[2], length.out = 100)
  data.frame(x = xs, y = predict(fit, newdata = data.frame(delta_age_old_minus_young = xs)))
}

global_range <- range(df_all$delta_age_old_minus_young)

# --- Combined overlay: solid over each clock's own data range, dashed
#     extrapolated out to the shared axis range for visual comparability ---
overlay_lines <- do.call(rbind, lapply(split(df_all, df_all$clock), function(sub) {
  own_range <- range(sub$delta_age_old_minus_young)
  pieces <- list(transform(fit_line(sub, own_range), clock = sub$clock[1], segment = "solid"))
  if (own_range[1] > global_range[1])
    pieces <- c(pieces, list(transform(fit_line(sub, c(global_range[1], own_range[1])), clock = sub$clock[1], segment = "dashed")))
  if (own_range[2] < global_range[2])
    pieces <- c(pieces, list(transform(fit_line(sub, c(own_range[2], global_range[2])), clock = sub$clock[1], segment = "dashed")))
  do.call(rbind, pieces)
}))

legend_labels <- setNames(
  sprintf("%s: slope=%.2f, r=%.2f, p=%.1e, n=%d", fit_stats$clock, fit_stats$slope, fit_stats$r, fit_stats$p, fit_stats$n),
  fit_stats$clock
)

p_overlay <- ggplot(df_all, aes(delta_age_old_minus_young, delta_reconst_minus_young, color = clock)) +
  geom_point(alpha = 0.55, size = 1.3) +
  geom_line(data = overlay_lines, aes(x, y, color = clock, linetype = segment), linewidth = 1.1) +
  scale_linetype_manual(values = c(solid = "solid", dashed = "dashed"), guide = "none") +
  scale_color_manual(values = clock_colors, labels = legend_labels, name = NULL) +
  geom_hline(yintercept = 0, color = "grey85") +
  geom_vline(xintercept = 0, color = "grey85") +
  labs(x = "Δ methylation: Old baseline − Young baseline (% pts, per CpG)",
       y = "Δ methylation: Young 10-cell reconst. − Young baseline (% pts, per CpG)",
       title = "Per-CpG change: natural aging vs. forced-proliferation transplant",
       subtitle = "GSE44117; dashed = extrapolated beyond that clock's own covered-CpG range") +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        panel.grid = element_blank())

ggsave(file.path(out_dir, "overlay_all_clocks.png"), p_overlay, width = 8.5, height = 7.5, dpi = 150)

# --- One separate plot per clock ---
for (cn in names(clock_colors)) {
  sub <- df_all[df_all$clock == cn, ]
  st  <- fit_stats[fit_stats$clock == cn, ]
  line_df <- fit_line(sub, range(sub$delta_age_old_minus_young))

  p <- ggplot(sub, aes(delta_age_old_minus_young, delta_reconst_minus_young)) +
    geom_point(color = clock_colors[cn], alpha = 0.6, size = 1.5) +
    geom_line(data = line_df, aes(x, y), color = clock_colors[cn], linewidth = 1.3) +
    geom_hline(yintercept = 0, color = "grey85") +
    geom_vline(xintercept = 0, color = "grey85") +
    annotate("label", x = -Inf, y = Inf, hjust = 0, vjust = 1,
             label = sprintf("slope=%.2f\nr=%.2f, p=%.1e\nn=%d", st$slope, st$r, st$p, st$n), size = 3.2) +
    labs(x = "Δ methylation: Old baseline − Young baseline (% pts, per CpG)",
         y = "Δ methylation: Young 10-cell reconst. − Young baseline (% pts, per CpG)",
         title = paste0(cn, " clock CpGs"),
         subtitle = "Natural aging vs. forced-proliferation transplant (GSE44117)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  ggsave(file.path(out_dir, paste0(tolower(cn), "_clock_correlation.png")), p, width = 6.5, height = 6, dpi = 150)
}

message("Done. All outputs written to: ", out_dir)
