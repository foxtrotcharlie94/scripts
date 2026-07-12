library(ggplot2)
library(dplyr)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

# ── 1. Load counts ────────────────────────────────────────────────────────────
counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

# ── 2. Compute ratios per cell ────────────────────────────────────────────────
compute_ratios <- function(counts, cell_type_label) {
  total_umis <- colSums(counts)
  sample_ids <- sub("_.*", "", colnames(counts))
  exposure   <- sub("[0-9]+$", "", sample_ids)
  
  trp53_gene <- grep("^Trp53$", rownames(counts), value = TRUE, ignore.case = TRUE)
  hsp_genes  <- grep("^Hsp",    rownames(counts), value = TRUE, ignore.case = TRUE)
  mt_genes   <- grep("^mt-",    rownames(counts), value = TRUE, ignore.case = TRUE)
  
  cat(cell_type_label, "- Trp53:", trp53_gene, "\n")
  cat(cell_type_label, "- HSP genes:", length(hsp_genes), "\n")
  cat(cell_type_label, "- MT genes:", length(mt_genes), "\n")
  
  data.frame(
    sample      = sample_ids,
    exposure    = exposure,
    cell_type   = cell_type_label,
    trp53_ratio = as.numeric(counts[trp53_gene, ]) / total_umis,
    hsp_ratio   = colSums(counts[hsp_genes, ])     / total_umis,
    mt_ratio    = colSums(counts[mt_genes, ])      / total_umis
  )
}

df_hsc  <- compute_ratios(counts_hsc,  "LT-HSC")
df_mono <- compute_ratios(counts_mono, "Classical Monocytes")
df_all  <- bind_rows(df_hsc, df_mono) %>%
  mutate(cell_type = factor(cell_type, levels = c("LT-HSC", "Classical Monocytes")))

# ── 3. Correlation helper ─────────────────────────────────────────────────────
cor_label <- function(x, y) {
  # Remove zeros for log-scale correlations
  keep  <- x > 0 & y > 0
  x_f   <- x[keep]
  y_f   <- y[keep]
  r     <- cor(log10(x_f), log10(y_f), method = "spearman")
  p     <- cor.test(log10(x_f), log10(y_f), method = "spearman")$p.value
  p_fmt <- ifelse(p < 0.001, "p<0.001", paste0("p=", signif(p, 2)))
  paste0("Spearman r=", round(r, 2), "\n", p_fmt,
         "\n(n=", sum(keep), " cells with >0 in both)")
}

# Compute correlation labels per cell type
cor_hsp_trp53_hsc  <- cor_label(df_hsc$hsp_ratio,  df_hsc$trp53_ratio)
cor_hsp_trp53_mono <- cor_label(df_mono$hsp_ratio, df_mono$trp53_ratio)
cor_hsp_mt_hsc     <- cor_label(df_hsc$hsp_ratio,  df_hsc$mt_ratio)
cor_hsp_mt_mono    <- cor_label(df_mono$hsp_ratio, df_mono$mt_ratio)

cor_labels_trp53 <- data.frame(
  cell_type = factor(c("LT-HSC", "Classical Monocytes"),
                     levels = c("LT-HSC", "Classical Monocytes")),
  label = c(cor_hsp_trp53_hsc, cor_hsp_trp53_mono)
)

cor_labels_mt <- data.frame(
  cell_type = factor(c("LT-HSC", "Classical Monocytes"),
                     levels = c("LT-HSC", "Classical Monocytes")),
  label = c(cor_hsp_mt_hsc, cor_hsp_mt_mono)
)

# ── 4. Plot: HSP ratio vs Trp53 ratio ─────────────────────────────────────────
p_hsp_trp53 <- ggplot(df_all %>% filter(hsp_ratio > 0, trp53_ratio > 0),
                      aes(x = hsp_ratio, y = trp53_ratio, color = exposure)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              aes(group = 1), color = "black") +
  geom_text(data = cor_labels_trp53,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3,
            size = 3, color = "grey20", fontface = "bold") +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  facet_wrap(~ cell_type, scales = "free") +
  labs(title    = "Heat shock vs Trp53 expression per cell",
       subtitle = "Each dot = one cell | X and Y axes log scale | Only cells with >0 UMIs in both genes shown",
       x = "HSP family UMIs / total UMIs",
       y = "Trp53 UMIs / total UMIs",
       color = "Exposure") +
  theme_classic() +
  theme(legend.position = "top",
        strip.text      = element_text(face = "bold", size = 11),
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 9, color = "grey40"))

# ── 5. Plot: HSP ratio vs MT ratio ────────────────────────────────────────────
p_hsp_mt <- ggplot(df_all %>% filter(hsp_ratio > 0, mt_ratio > 0),
                   aes(x = hsp_ratio, y = mt_ratio, color = exposure)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              aes(group = 1), color = "black") +
  geom_text(data = cor_labels_mt,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3,
            size = 3, color = "grey20", fontface = "bold") +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  facet_wrap(~ cell_type, scales = "free") +
  labs(title    = "Heat shock vs Mitochondrial gene expression per cell",
       subtitle = "Each dot = one cell | X and Y axes log scale | Only cells with >0 UMIs in both shown",
       x = "HSP family UMIs / total UMIs",
       y = "Mitochondrial UMIs / total UMIs",
       color = "Exposure") +
  theme_classic() +
  theme(legend.position = "top",
        strip.text      = element_text(face = "bold", size = 11),
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 9, color = "grey40"))

# ── 6. Combined ───────────────────────────────────────────────────────────────
p_combined <- p_hsp_trp53 / p_hsp_mt

ggsave("hsp_vs_trp53_per_cell.pdf",    plot = p_hsp_trp53, width = 12, height = 6)
ggsave("hsp_vs_mt_per_cell.pdf",       plot = p_hsp_mt,    width = 12, height = 6)
ggsave("hsp_correlations_combined.pdf",plot = p_combined,  width = 12, height = 12)

cat("Correlation plots saved!\n")
