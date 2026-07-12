## ============================================================
##  HSC stochastic division model — all plots
##  Outputs:
##    HSC_panel.pdf / .png          — 8-panel combined figure
##    HSC_01_division_burden.pdf/png
##    HSC_02_pool_size.pdf/png
##    HSC_03_FRAC_division_burden.pdf/png
##    HSC_04_A2_division_burden.pdf/png
##    HSC_05_FRAC_pool_size.pdf/png
##    HSC_06_A2_division_burden_comparison.pdf/png
##    HSC_07_burden_distribution.pdf/png
##    HSC_08_A2_pool_size_logscale.pdf/png
## ============================================================

library(tidyverse)
library(patchwork)
library(scales)

set.seed(42)

## ── Parameters ────────────────────────────────────────────────────────────────
N0        <- 1000
TIMESTEPS <- 100
FRAC      <- 0.25
A2        <- 0.15
A1        <- 0.75
B         <- 0.10
NSIM      <- 50
NSIM_HIST <- 500

stopifnot("A2+A1+B must equal 1" = abs(A2+A1+B-1) < 1e-9)

a2_eff <- FRAC * A2
a1_eff <- FRAC * A1
b_eff  <- FRAC * B

## ── Simulation function ───────────────────────────────────────────────────────
run_sim <- function(n0, timesteps, frac, A2, A1, B) {
  division_numbers   <- c(n0, rep(0L, timesteps))
  pool_size          <- integer(timesteps + 1)
  div_burden_mean    <- numeric(timesteps + 1)
  pool_size[1]       <- n0
  div_burden_mean[1] <- 0
  
  for (i in seq_len(timesteps)) {
    total_pool <- sum(division_numbers)
    if (total_pool == 0) {
      pool_size[i+1]       <- 0L
      div_burden_mean[i+1] <- div_burden_mean[i]
      next
    }
    n_dividing <- rbinom(n=timesteps+1, size=division_numbers, prob=frac)
    selfrenew  <- rbinom(n=timesteps+1, size=n_dividing,  prob=A2)
    remaining  <- n_dividing - selfrenew
    asym       <- rbinom(n=timesteps+1, size=remaining, prob=A1/(A1+B))
    death      <- remaining - asym
    
    division_numbers <- division_numbers - n_dividing +
      c(0L, n_dividing[-(timesteps+1)]) +
      c(0L, selfrenew[-(timesteps+1)]) -
      c(0L, death[-(timesteps+1)])
    division_numbers <- pmax(0L, division_numbers)
    
    total_new <- sum(division_numbers)
    pool_size[i+1] <- total_new
    div_burden_mean[i+1] <- if (total_new > 0)
      weighted.mean(0:timesteps, division_numbers) else div_burden_mean[i]
  }
  list(pool=pool_size, div_burden=div_burden_mean)
}

## ── Run simulations ───────────────────────────────────────────────────────────
message("Running ", NSIM, " simulations...")
sims <- lapply(seq_len(NSIM), function(i)
  run_sim(N0, TIMESTEPS, FRAC, A2, A1, B))

message("Running ", NSIM_HIST, " simulations for histogram...")
sims_hist <- lapply(seq_len(NSIM_HIST), function(i)
  run_sim(N0, TIMESTEPS, FRAC, A2, A1, B))

t_seq      <- 0:TIMESTEPS
analytic_b <- (2*a2_eff+a1_eff)*t_seq/(1+a2_eff-b_eff)
analytic_p <- N0*(1+a2_eff-b_eff)^t_seq

## ── Trajectory data frames ────────────────────────────────────────────────────
burden_df <- do.call(rbind, lapply(seq_along(sims), function(i)
  data.frame(sim=i, t=t_seq, burden=sims[[i]]$div_burden)))
pool_df   <- do.call(rbind, lapply(seq_along(sims), function(i)
  data.frame(sim=i, t=t_seq, pool=sims[[i]]$pool)))

burden_sum <- burden_df %>% group_by(t) %>%
  summarise(mean=mean(burden), q05=quantile(burden,.05), q95=quantile(burden,.95), .groups="drop")
pool_sum   <- pool_df %>% group_by(t) %>%
  summarise(mean=mean(pool), q05=quantile(pool,.05), q95=quantile(pool,.95), .groups="drop")

analytic_df <- data.frame(t=t_seq, burden=analytic_b, pool=analytic_p)

## ── Colour helpers ────────────────────────────────────────────────────────────
# Generate n colours from a palette avoiding near-white ends
pal <- function(palette, n, start=0.15, end=0.95)
  colorRampPalette(RColorBrewer::brewer.pal(max(3,min(9,n)), palette))(n)

# Safe wrapper if RColorBrewer not available — use viridis
make_pal <- function(palette, n) {
  tryCatch(pal(palette, n), error=function(e)
    viridis::viridis(n, option=switch(palette,
                                      Blues="D", Greens="E", Oranges="A", "D")))
}

frac_vals <- seq(0.05, 0.50, by=0.05)
A2_vals   <- seq(0.10, 0.50, by=0.05)

## ── Theme ─────────────────────────────────────────────────────────────────────
th <- theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, color="grey40"),
        legend.position="right",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

subtitle_main <- sprintf("frac=%.2f  A2=%.2f  A1=%.2f  B=%.2f  N0=%d  n=%d sims  |  shaded=5–95th pct",
                         FRAC, A2, A1, B, N0, NSIM)

## ── Plot 1: Division burden sim vs analytical ─────────────────────────────────
p1 <- ggplot(burden_sum, aes(x=t)) +
  geom_ribbon(aes(ymin=q05, ymax=q95), fill="#534AB7", alpha=0.18) +
  geom_line(aes(y=mean, color="simulation (mean)"), lwd=1.2) +
  geom_line(data=analytic_df, aes(x=t, y=burden, color="analytical"),
            lwd=1.1, lty="dashed") +
  scale_color_manual(values=c("simulation (mean)"="#534AB7","analytical"="#C0392B"), name=NULL) +
  labs(title="Mean division burden per HSC", subtitle=subtitle_main,
       x="Timestep", y="Mean divisions per cell") + th

## ── Plot 2: Pool size sim vs analytical ───────────────────────────────────────
p2 <- ggplot(pool_sum, aes(x=t)) +
  geom_ribbon(aes(ymin=q05, ymax=q95), fill="#378ADD", alpha=0.18) +
  geom_line(aes(y=mean, color="simulation (mean)"), lwd=1.2) +
  geom_line(data=analytic_df, aes(x=t, y=pool, color="N0·(1+a2_eff−b_eff)^t"),
            lwd=1.1, lty="dashed") +
  geom_hline(yintercept=N0, color="grey60", lty="dotted") +
  scale_color_manual(values=c("simulation (mean)"="#378ADD",
                              "N0·(1+a2_eff−b_eff)^t"="#C0392B"), name=NULL) +
  scale_y_continuous(labels=comma) +
  labs(title="HSC pool size N(t)", x="Timestep", y="Pool size") + th

## ── Plot 3: FRAC → division burden ───────────────────────────────────────────
frac_burden_df <- do.call(rbind, lapply(frac_vals, function(fv) {
  ae2=fv*A2; ae1=fv*A1; be=fv*B
  data.frame(t=t_seq, frac=fv,
             burden=(2*ae2+ae1)*t_seq/(1+ae2-be))
}))

p3 <- ggplot(frac_burden_df, aes(x=t, y=burden, group=factor(frac), color=frac)) +
  geom_line(lwd=0.9) +
  scale_color_distiller(palette="Blues", direction=1, name="FRAC",
                        limits=c(0.05,0.50)) +
  labs(title="Impact of FRAC on division burden",
       subtitle=sprintf("A2=%.2f  A1=%.2f  B=%.2f fixed", A2, A1, B),
       x="Timestep", y="Mean division burden (analytical)") + th

## ── Plot 4: A2 → division burden ──────────────────────────────────────────────
A2_burden_df <- do.call(rbind, lapply(A2_vals, function(a2v) {
  ae2=FRAC*a2v; ae1=FRAC*max(0,1-a2v-B); be=FRAC*B
  data.frame(t=t_seq, A2=a2v,
             burden=(2*ae2+ae1)*t_seq/(1+ae2-be))
}))

p4 <- ggplot(A2_burden_df, aes(x=t, y=burden, group=factor(A2), color=A2)) +
  geom_line(lwd=0.9) +
  scale_color_distiller(palette="Greens", direction=1, name="A2",
                        limits=c(0.10,0.50)) +
  labs(title="Impact of A2 on division burden",
       subtitle=sprintf("frac=%.2f  B=%.2f fixed  |  A1 = 1 − A2 − B", FRAC, B),
       x="Timestep", y="Mean division burden (analytical)") + th

## ── Plot 5: FRAC → pool size (highlighted) ────────────────────────────────────
frac_pool_df <- do.call(rbind, lapply(frac_vals, function(fv)
  data.frame(t=t_seq, frac=fv, pool=N0*(1+fv*A2-fv*B)^t_seq)))

highlight_df <- data.frame(
  t=rep(t_seq,2),
  pool=c(N0*(1+0.25*A2-0.25*B)^t_seq, N0*(1+0.50*A2-0.50*B)^t_seq),
  grp=rep(c("frac=0.25","frac=0.50"),each=length(t_seq))
)

p5 <- ggplot(frac_pool_df, aes(x=t, y=pool, group=factor(frac), color=frac)) +
  geom_line(lwd=0.9) +
  geom_line(data=filter(highlight_df, grp=="frac=0.25"),
            aes(x=t, y=pool, group=1, linetype="frac=0.25"),
            color="#C0392B", lwd=2, inherit.aes=FALSE) +
  geom_line(data=filter(highlight_df, grp=="frac=0.50"),
            aes(x=t, y=pool, group=1, linetype="frac=0.50"),
            color="#7B2FBE", lwd=2, inherit.aes=FALSE) +
  geom_hline(yintercept=N0, color="grey60", lty="dotted") +
  scale_color_distiller(palette="Blues", direction=1, name="FRAC",
                        limits=c(0.05,0.50)) +
  scale_linetype_manual(
    name   = "Highlighted",
    values = c("frac=0.25"="solid","frac=0.50"="solid"),
    guide  = guide_legend(
      override.aes=list(color=c("#C0392B","#7B2FBE"), lwd=c(2,2))
    )
  ) +
  scale_y_continuous(labels=comma) +
  labs(title="Impact of FRAC on HSC pool size",
       subtitle=sprintf("A2=%.2f  B=%.2f fixed  |  FRAC range 0.05–0.50", A2, B),
       x="Timestep", y="Pool size (analytical)") + th

## ── Plot 6: A2 → division burden frac=0.25 vs 0.50 ───────────────────────────
A2_cmp_df <- do.call(rbind, lapply(A2_vals, function(a2v) {
  rbind(
    data.frame(t=t_seq, A2=a2v, frac_grp="frac=0.25",
               burden=(2*0.25*a2v+0.25*max(0,1-a2v-B))*t_seq/
                 (1+0.25*a2v-0.25*B)),
    data.frame(t=t_seq, A2=a2v, frac_grp="frac=0.50",
               burden=(2*0.50*a2v+0.50*max(0,1-a2v-B))*t_seq/
                 (1+0.50*a2v-0.50*B))
  )
}))

p6 <- ggplot(A2_cmp_df, aes(x=t, y=burden, group=interaction(A2, frac_grp),
                            color=A2, lty=frac_grp)) +
  geom_line(lwd=0.9, alpha=0.85) +
  scale_color_distiller(palette="Greens", direction=1, name="A2",
                        limits=c(0.10,0.50)) +
  scale_linetype_manual(values=c("frac=0.25"="solid","frac=0.50"="dashed"),
                        name="FRAC") +
  labs(title="A2 on division burden: frac=0.25 vs frac=0.50",
       subtitle=sprintf("B=%.2f fixed  |  A1=1−A2−B  |  A2 range 0.10–0.50", B),
       x="Timestep", y="Mean division burden (analytical)") + th

## ── Plot 7: Division burden distribution (n=500) ──────────────────────────────
final_burden <- sapply(sims_hist, function(s) s$div_burden[TIMESTEPS+1])
dist_df <- data.frame(burden=final_burden)

p7 <- ggplot(dist_df, aes(x=burden)) +
  geom_histogram(bins=40, fill="#1D9E75", color="white", lwd=0.3) +
  geom_vline(xintercept=mean(final_burden), color="#534AB7", lty="dashed", lwd=1.2) +
  geom_vline(xintercept=analytic_b[TIMESTEPS+1], color="#C0392B", lty="dashed", lwd=1.2) +
  annotate("text", x=mean(final_burden), y=Inf, vjust=2, hjust=-0.1,
           label=sprintf("sim mean\n%.2f", mean(final_burden)),
           color="#534AB7", size=3.5) +
  annotate("text", x=analytic_b[TIMESTEPS+1], y=Inf, vjust=2, hjust=1.1,
           label=sprintf("analytical\n%.2f", analytic_b[TIMESTEPS+1]),
           color="#C0392B", size=3.5) +
  labs(title=sprintf("Division burden distribution at t=%d  (n=%d)", TIMESTEPS, NSIM_HIST),
       x="Mean divisions per HSC", y="Simulations") + th

## ── Plot 8: A2 → pool size frac=0.25 vs 0.50, log scale ──────────────────────
A2_pool_cmp_df <- do.call(rbind, lapply(A2_vals, function(a2v) {
  rbind(
    data.frame(t=t_seq, A2=a2v, frac_grp="frac=0.25",
               pool=N0*(1+0.25*a2v-0.25*B)^t_seq),
    data.frame(t=t_seq, A2=a2v, frac_grp="frac=0.50",
               pool=N0*(1+0.50*a2v-0.50*B)^t_seq)
  )
}))

# Use two separate colour scales via ggnewscale if available, else facet
p8 <- ggplot(A2_pool_cmp_df,
             aes(x=t, y=pool, group=interaction(A2, frac_grp),
                 color=A2, lty=frac_grp)) +
  geom_line(lwd=0.9) +
  geom_hline(yintercept=N0, color="grey60", lty="dotted") +
  scale_color_distiller(palette="Oranges", direction=1, name="A2",
                        limits=c(0.10,0.50)) +
  scale_linetype_manual(values=c("frac=0.25"="solid","frac=0.50"="dashed"),
                        name="FRAC") +
  scale_y_log10(labels=comma) +
  labs(title="A2 on pool size: frac=0.25 (solid) vs frac=0.50 (dashed)",
       subtitle=sprintf("B=%.2f fixed  |  A2 range 0.10–0.50  |  log y-axis", B),
       x="Timestep", y="Pool size — log scale (analytical)") + th

## ── Save individual plots ─────────────────────────────────────────────────────
plot_list <- list(
  "HSC_01_division_burden"        = p1,
  "HSC_02_pool_size"              = p2,
  "HSC_03_FRAC_division_burden"   = p3,
  "HSC_04_A2_division_burden"     = p4,
  "HSC_05_FRAC_pool_size"         = p5,
  "HSC_06_A2_burden_comparison"   = p6,
  "HSC_07_burden_distribution"    = p7,
  "HSC_08_A2_pool_size_logscale"  = p8
)

for (nm in names(plot_list)) {
  ggsave(paste0(nm, ".pdf"), plot_list[[nm]], width=8, height=5.5, limitsize=FALSE)
  ggsave(paste0(nm, ".png"), plot_list[[nm]], width=8, height=5.5, dpi=150, limitsize=FALSE)
  message("Saved: ", nm)
}

## ── Combined 8-panel figure ───────────────────────────────────────────────────
message("\nBuilding combined panel...")

panel <- (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | p8) +
  plot_annotation(
    title    = sprintf("HSC stochastic division model  |  frac=%.2f  A2=%.2f  A1=%.2f  B=%.2f  N0=%d",
                       FRAC, A2, A1, B, N0),
    theme    = theme(plot.title=element_text(face="bold", size=14, hjust=0.5))
  )

ggsave("HSC_panel.pdf", panel, width=16, height=22, limitsize=FALSE)
ggsave("HSC_panel.png", panel, width=16, height=22, dpi=120, limitsize=FALSE)

message("\nDone. Files written:")
message("  HSC_panel.pdf / .png")
for (nm in names(plot_list)) message("  ", nm, ".pdf / .png")

