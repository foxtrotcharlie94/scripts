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

