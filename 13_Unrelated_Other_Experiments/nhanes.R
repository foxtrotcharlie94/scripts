# ============================================================
# NHANES CBC: Download & plot all parameters vs age (M/F)
# Pools multiple cycles; uses survey-weighted LOESS smoothing
# ============================================================

# ── 0. Install / load packages ────────────────────────────────
pkgs <- c("nhanesA", "dplyr", "tidyr", "ggplot2", "patchwork", "survey", "quantreg")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(nhanesA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(survey)
library(quantreg)

# Explicitly recover dplyr generics masked by survey/MASS
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
summarise <- dplyr::summarise

# ── 1. Define cycles to pull ─────────────────────────────────
# Using 2005–2023 (post-instrument-standardisation era).
# Exclude 2013-2014 (CBC_H) trend artefacts if you want;
# keep here for completeness but flag in notes.
cycles <- list(
  list(cbc = "CBC_D", demo = "DEMO_D", yr = "2005-06"),
  list(cbc = "CBC_E", demo = "DEMO_E", yr = "2007-08"),
  list(cbc = "CBC_F", demo = "DEMO_F", yr = "2009-10"),
  list(cbc = "CBC_G", demo = "DEMO_G", yr = "2011-12"),
  list(cbc = "CBC_H", demo = "DEMO_H", yr = "2013-14"),  # instrument change note
  list(cbc = "CBC_I", demo = "DEMO_I", yr = "2015-16"),
  list(cbc = "P_CBC",  demo = "P_DEMO",  yr = "2017-Mar2020"),  # COVID truncation: 2017-18 + partial 2019-20 combined by CDC
  list(cbc = "CBC_L", demo = "DEMO_L", yr = "2021-23")
)

# ── 2. CBC variables of interest ─────────────────────────────
cbc_vars <- c(
  "LBXWBCSI",  # WBC total (10^3/uL)
  "LBXLYPCT",  # Lymphocyte %
  "LBDLYMNO",  # Lymphocyte absolute (10^3/uL)
  "LBXMOPCT",  # Monocyte %
  "LBDMONO",   # Monocyte absolute (10^3/uL)
  "LBXNEPCT",  # Neutrophil %
  "LBDNENO",   # Neutrophil absolute (10^3/uL)
  "LBXEOPCT",  # Eosinophil %
  "LBDEONO",   # Eosinophil absolute (10^3/uL)
  "LBXBAPCT",  # Basophil %
  "LBDBANO",   # Basophil absolute (10^3/uL)
  "LBXHGB",    # Hemoglobin (g/dL)
  "LBXHCT",    # Hematocrit (%)
  "LBXRBCSI",  # RBC (10^6/uL)
  "LBXMCVSI",  # MCV (fL)
  "LBXMCHSI",  # MCH (pg)
  "LBXMC",     # MCHC (g/dL)  -- note: LBXMCHSI in some cycles
  "LBXRDW",    # RDW (%)
  "LBXPLTSI",  # Platelets (10^3/uL)
  "LBXMPSI"    # MPV (fL)
)

demo_vars <- c("SEQN", "RIDAGEYR", "RIAGENDR", "WTMEC2YR",
               "SDMVPSU", "SDMVSTRA", "RIDRETH1")

# ── 3. Download and merge ─────────────────────────────────────
message("Downloading NHANES data (this may take a few minutes)...")

all_data <- lapply(cycles, function(cyc) {
  message("  Fetching cycle: ", cyc$yr)
  
  cbc_raw  <- tryCatch(nhanes(cyc$cbc,  translated = FALSE), error = function(e) NULL)
  demo_raw <- tryCatch(nhanes(cyc$demo, translated = FALSE), error = function(e) NULL)
  
  if (is.null(cbc_raw) || is.null(demo_raw)) {
    message("    [SKIPPED - download failed]")
    return(NULL)
  }
  
  # Keep only available vars
  cbc_keep  <- intersect(c("SEQN", cbc_vars), names(cbc_raw))
  demo_keep <- intersect(demo_vars, names(demo_raw))
  
  merged <- inner_join(
    cbc_raw[, cbc_keep],
    demo_raw[, demo_keep],
    by = "SEQN"
  ) %>%
    mutate(cycle = cyc$yr)
  
  merged
})

dat <- bind_rows(Filter(Negate(is.null), all_data))

message("Total rows after merging: ", nrow(dat))
message("Columns in dat: ", paste(sort(names(dat)), collapse = ", "))

# Defensive check — rename if nhanesA returned translated column names
if (!"RIDAGEYR" %in% names(dat) && "Age at Screening Adjudicated - Recode" %in% names(dat))
  dat <- rename(dat, RIDAGEYR = `Age at Screening Adjudicated - Recode`)
if (!"RIAGENDR" %in% names(dat) && "Gender" %in% names(dat))
  dat <- rename(dat, RIAGENDR = Gender)
if (!"WTMEC2YR" %in% names(dat) && "Full Sample 2 Year MEC Exam Weight" %in% names(dat))
  dat <- rename(dat, WTMEC2YR = `Full Sample 2 Year MEC Exam Weight`)

# Fallback: if still missing, try fetching demo with nhanesTranslate = FALSE
if (!"RIDAGEYR" %in% names(dat)) {
  message("RIDAGEYR still missing — re-fetching demo files untranslated...")
  all_data2 <- lapply(cycles, function(cyc) {
    demo_raw <- tryCatch(nhanes(cyc$demo, translated = FALSE), error = function(e) NULL)
    if (is.null(demo_raw)) return(NULL)
    demo_keep <- intersect(demo_vars, names(demo_raw))
    demo_raw[, demo_keep] %>% mutate(cycle = cyc$yr)
  })
  demo_all <- bind_rows(Filter(Negate(is.null), all_data2))
  # Re-merge with CBC columns already in dat
  cbc_cols <- setdiff(names(dat), demo_vars)
  cbc_cols <- c("SEQN", intersect(cbc_cols, names(dat)))
  dat <- inner_join(dat[, intersect(cbc_cols, names(dat))],
                    demo_all, by = "SEQN")
  message("Columns after re-fetch: ", paste(sort(names(dat)), collapse = ", "))
}

# ── 4. Clean / recode ─────────────────────────────────────────
# nhanesA translates coded values by default; handle both string and numeric
message("Unique RIAGENDR values: ", paste(unique(dat$RIAGENDR), collapse = ", "))

dat <- dat %>%
  filter(RIDAGEYR >= 6, WTMEC2YR > 0) %>%
  mutate(
    sex = case_when(
      RIAGENDR %in% c(1, "1", "Male")   ~ "Male",
      RIAGENDR %in% c(2, "2", "Female") ~ "Female",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Male", "Female")),
    age = as.numeric(RIDAGEYR)
  ) %>%
  filter(!is.na(sex), !is.na(age))

message("Rows after sex/age filter: ", nrow(dat))

# Handle MCH/MCHC column name variation across cycles
if ("LBXMCH" %in% names(dat) && !"LBXMCHSI" %in% names(dat)) {
  dat <- rename(dat, LBXMCHSI = LBXMCH)
} else if ("LBXMCH" %in% names(dat) && "LBXMCHSI" %in% names(dat)) {
  dat$LBXMCHSI <- coalesce(dat$LBXMCHSI, dat$LBXMCH)
}
if ("LBXMC" %in% names(dat)) {
  if (!"LBXMCSI" %in% names(dat)) {
    dat <- rename(dat, LBXMCSI = LBXMC)
  } else {
    dat$LBXMCSI <- coalesce(dat$LBXMCSI, dat$LBXMC)
  }
}

# ── 4b. Comprehensive exclusion of conditions affecting CBC ───
#
# Table suffix lookup: CBC_D -> suffix "D", P_CBC -> suffix "P" etc.
suffix_map <- list(
  "CBC_D" = "D", "CBC_E" = "E", "CBC_F" = "F", "CBC_G" = "G",
  "CBC_H" = "H", "CBC_I" = "I", "P_CBC"  = "P", "CBC_L" = "L"
)

safe_nhanes <- function(tbl) {
  tryCatch(nhanes(tbl, translated = FALSE),
           error = function(e) { message("  [missing] ", tbl); NULL })
}

flag_yes <- function(vec) vec %in% c(1, "1", "Yes", "YES")

`%||%` <- function(a, b) if (!is.null(a)) a else b

excl_seqn <- c()
n_excl    <- list()  # track counts per reason for reporting

for (cyc in cycles) {
  sfx <- suffix_map[[cyc$cbc]]
  message("  Exclusion downloads for cycle: ", cyc$yr)
  
  # ── MCQ: Medical Conditions ──────────────────────────────────
  # MCQ220 = cancer (any)
  # MCQ092 = told have anaemia
  # MCQ160L = liver condition
  # MCQ160A = arthritis (proxy for autoimmune — imperfect but standard)
  # MCQ160N = lupus
  mcq <- safe_nhanes(if (sfx == "P") "P_MCQ" else paste0("MCQ_", sfx))
  if (!is.null(mcq)) {
    for (v in c("MCQ220", "MCQ092", "MCQ160L", "MCQ160A", "MCQ160N")) {
      if (v %in% names(mcq)) {
        bad <- mcq$SEQN[flag_yes(mcq[[v]])]
        excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
        n_excl[[v]] <- (n_excl[[v]] %||% 0) + sum(!is.na(bad))
      }
    }
  }
  
  # ── KIQ: Kidney Conditions ───────────────────────────────────
  # KIQ022 = told have weak/failing kidneys
  # KIQ025 = currently on dialysis
  # File named KIQ_U_ in most cycles, KIQ_ in a few early ones
  kiq_name <- if (sfx == "P") "P_KIQ_U" else paste0("KIQ_U_", sfx)
  kiq <- safe_nhanes(kiq_name)
  if (is.null(kiq))
    kiq <- safe_nhanes(if (sfx == "P") "P_KIQ" else paste0("KIQ_", sfx))
  if (!is.null(kiq)) {
    for (v in c("KIQ022", "KIQ025")) {
      if (v %in% names(kiq)) {
        bad <- kiq$SEQN[flag_yes(kiq[[v]])]
        excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
        n_excl[[v]] <- (n_excl[[v]] %||% 0) + sum(!is.na(bad))
      }
    }
  }
  
  # ── RHQ: Reproductive Health (pregnancy) ─────────────────────
  # RHD143 = currently pregnant (confirmed)
  # RHQ141 = thinks might be pregnant
  rhq <- safe_nhanes(if (sfx == "P") "P_RHQ" else paste0("RHQ_", sfx))
  if (!is.null(rhq)) {
    for (v in c("RHD143", "RHQ141")) {
      if (v %in% names(rhq)) {
        bad <- rhq$SEQN[flag_yes(rhq[[v]])]
        excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
        n_excl[[v]] <- (n_excl[[v]] %||% 0) + sum(!is.na(bad))
      }
    }
  }
  
  # ── HIV (lab) ────────────────────────────────────────────────
  # LBXHIV = HIV antibody result (1 = reactive/positive)
  hiv <- safe_nhanes(if (sfx == "P") "P_HIV" else paste0("HIV_", sfx))
  if (!is.null(hiv) && "LBXHIV" %in% names(hiv)) {
    bad <- hiv$SEQN[hiv$LBXHIV %in% c(1, "Reactive", "Positive")]
    excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
    n_excl[["LBXHIV"]] <- (n_excl[["LBXHIV"]] %||% 0) + sum(!is.na(bad))
  }
  
  # ── CRP: Acute inflammation proxy ────────────────────────────
  # LBXCRP > 10 mg/L = likely acute infection/inflammation
  # File named CRP_ for 2005-2010 cycles, HSCRP_ from 2011 onwards
  crp_name <- if (sfx %in% c("D", "E", "F")) paste0("CRP_", sfx) else
    if (sfx == "P") "P_HSCRP" else paste0("HSCRP_", sfx)
  crp <- safe_nhanes(crp_name)
  if (!is.null(crp)) {
    crp_var <- intersect(c("LBXHSCRP", "LBXCRP"), names(crp))[1]
    if (!is.na(crp_var)) {
      bad <- crp$SEQN[!is.na(crp[[crp_var]]) & crp[[crp_var]] > 10]
      excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
      n_excl[["CRP>10"]] <- (n_excl[["CRP>10"]] %||% 0) + sum(!is.na(bad))
    }
  }
  
  # ── SMQ: Smoking ─────────────────────────────────────────────
  # SMQ020 = smoked at least 100 cigarettes in life (ever-smoker)
  smq <- safe_nhanes(if (sfx == "P") "P_SMQ" else paste0("SMQ_", sfx))
  if (!is.null(smq) && "SMQ020" %in% names(smq)) {
    bad <- smq$SEQN[flag_yes(smq$SMQ020)]
    excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
    n_excl[["SMQ020"]] <- (n_excl[["SMQ020"]] %||% 0) + sum(!is.na(bad))
  }
  
  # ── BMX: Obesity (BMI >= 35) ─────────────────────────────────
  bmx <- safe_nhanes(if (sfx == "P") "P_BMX" else paste0("BMX_", sfx))
  if (!is.null(bmx) && "BMXBMI" %in% names(bmx)) {
    bad <- bmx$SEQN[!is.na(bmx$BMXBMI) & bmx$BMXBMI >= 35]
    excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
    n_excl[["BMI>=35"]] <- (n_excl[["BMI>=35"]] %||% 0) + sum(!is.na(bad))
  }
  
  # ── ALQ: Heavy alcohol (> 14 drinks/week) ────────────────────
  # ALQ130 = avg drinks per day in past 12 months -> *7 for weekly
  alq <- safe_nhanes(if (sfx == "P") "P_ALQ" else paste0("ALQ_", sfx))
  if (!is.null(alq) && "ALQ130" %in% names(alq)) {
    bad <- alq$SEQN[!is.na(alq$ALQ130) & alq$ALQ130 * 7 > 14]
    excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
    n_excl[["ALQ>14/wk"]] <- (n_excl[["ALQ>14/wk"]] %||% 0) + sum(!is.na(bad))
  }
  
  # ── DIQ: Diabetes ────────────────────────────────────────────
  # DIQ010: "told by doctor you have diabetes" (1 = Yes)
  diq <- safe_nhanes(if (sfx == "P") "P_DIQ" else paste0("DIQ_", sfx))
  if (!is.null(diq) && "DIQ010" %in% names(diq)) {
    bad <- diq$SEQN[flag_yes(diq$DIQ010)]
    excl_seqn <- c(excl_seqn, bad[!is.na(bad)])
    n_excl[["DIQ010"]] <- (n_excl[["DIQ010"]] %||% 0) + sum(!is.na(bad))
  }
}

excl_seqn <- unique(excl_seqn)
message("\n── Exclusion summary ──────────────────────────────")
for (nm in names(n_excl))
  message(sprintf("  %-15s : %d flagged across all cycles", nm, n_excl[[nm]]))
message(sprintf("  %-15s : %d unique SEQNs excluded", "TOTAL (unique)", length(excl_seqn)))
message("────────────────────────────────────────────────────\n")

dat <- dat %>% filter(!SEQN %in% excl_seqn)
message("Rows after exclusions: ", nrow(dat))

# ── 5. Define parameter metadata ─────────────────────────────
param_meta <- tribble(
  ~var,        ~label,                ~unit,
  "LBXWBCSI",  "WBC",                 "10³/µL",
  "LBXLYPCT",  "Lymphocyte %",        "%",
  "LBDLYMNO",  "Lymphocyte (abs)",    "10³/µL",
  "LBXMOPCT",  "Monocyte %",          "%",
  "LBDMONO",   "Monocyte (abs)",      "10³/µL",
  "LBXNEPCT",  "Neutrophil %",        "%",
  "LBDNENO",   "Neutrophil (abs)",    "10³/µL",
  "LBXEOPCT",  "Eosinophil %",        "%",
  "LBDEONO",   "Eosinophil (abs)",    "10³/µL",
  "LBXBAPCT",  "Basophil %",          "%",
  "LBDBANO",   "Basophil (abs)",      "10³/µL",
  "LBXHGB",    "Hemoglobin",          "g/dL",
  "LBXHCT",    "Hematocrit",          "%",
  "LBXRBCSI",  "RBC",                 "10⁶/µL",
  "LBXMCVSI",  "MCV",                 "fL",
  "LBXMCHSI",  "MCH",                 "pg",
  "LBXRDW",    "RDW",                 "%",
  "LBXPLTSI",  "Platelets",           "10³/µL",
  "LBXMPSI",   "MPV",                 "fL"
)
param_meta <- param_meta %>% filter(var %in% names(dat))
message("Parameters found: ", nrow(param_meta))

# ── 6. Reshape to long + per-variable outlier removal ─────────
dat_long <- dat %>%
  select(SEQN, age, sex, WTMEC2YR, all_of(param_meta$var)) %>%
  pivot_longer(cols = all_of(param_meta$var),
               names_to = "var", values_to = "value") %>%
  filter(!is.na(value)) %>%
  left_join(param_meta, by = "var") %>%
  # Remove extreme outliers per variable: outside Tukey 3×IQR fence
  group_by(var) %>%
  mutate(
    lo = quantile(value, 0.25, na.rm = TRUE) - 3 * IQR(value, na.rm = TRUE),
    hi = quantile(value, 0.75, na.rm = TRUE) + 3 * IQR(value, na.rm = TRUE)
  ) %>%
  filter(value >= lo, value <= hi) %>%
  select(-lo, -hi) %>%
  ungroup()

# ── 7. n per sex for labelling ────────────────────────────────
n_per_sex <- dat_long %>%
  filter(var == param_meta$var[1]) %>%
  count(sex) %>%
  mutate(lab = paste0(sex, ": n=", format(n, big.mark = ","))) %>%
  pull(lab) %>%
  paste(collapse = "  |  ")

message("Sample sizes: ", n_per_sex)

# ── 8. Plot function — all smooths, no bins ───────────────────
# Median line : weighted LOESS
# IQR ribbon  : quantile regression spline (rqss) at tau = 0.25 / 0.75
sex_colours <- c("Male" = "#2c7bb6", "Female" = "#d7191c")
sex_fills   <- c("Male" = "#abd9e9", "Female" = "#fdae61")

plot_param <- function(pm) {
  d <- dat_long %>% filter(var == pm$var)
  if (nrow(d) == 0) return(NULL)
  
  n_lab <- d %>%
    count(sex) %>%
    mutate(lab = paste0(sex, ": n=", format(n, big.mark = ","))) %>%
    pull(lab) %>%
    paste(collapse = "   ")
  
  # Compute smooth Q25/Q75 via rqss per sex, return prediction grid
  age_grid <- seq(6, 85, by = 1)
  ribbon_df <- lapply(levels(d$sex), function(s) {
    ds <- d %>% filter(sex == s)
    fit25 <- tryCatch(
      predict(rqss(value ~ qss(age, lambda = 5), tau = 0.25, data = ds),
              newdata = data.frame(age = age_grid)),
      error = function(e) rep(NA_real_, length(age_grid))
    )
    fit75 <- tryCatch(
      predict(rqss(value ~ qss(age, lambda = 5), tau = 0.75, data = ds),
              newdata = data.frame(age = age_grid)),
      error = function(e) rep(NA_real_, length(age_grid))
    )
    data.frame(age = age_grid, q25 = fit25, q75 = fit75,
               sex = factor(s, levels = levels(d$sex)))
  })
  ribbon_df <- bind_rows(ribbon_df)
  
  ggplot(d, aes(x = age, y = value, colour = sex, fill = sex)) +
    geom_ribbon(data = ribbon_df,
                aes(x = age, ymin = q25, ymax = q75, fill = sex, group = sex),
                inherit.aes = FALSE, alpha = 0.18) +
    geom_smooth(aes(weight = WTMEC2YR),
                method = "loess", formula = y ~ x,
                span = 0.4, se = FALSE, linewidth = 1.0) +
    scale_colour_manual(values = sex_colours, name = NULL) +
    scale_fill_manual(values = sex_fills, name = NULL) +
    scale_x_continuous(breaks = seq(10, 80, 10), limits = c(6, 85)) +
    labs(title = pm$label, subtitle = n_lab,
         x = "Age (years)", y = pm$unit) +
    theme_bw(base_size = 10) +
    theme(plot.title      = element_text(face = "bold", size = 10),
          plot.subtitle   = element_text(size = 7, colour = "grey40"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          panel.grid.minor = element_blank())
}

# ── 9. Build & save multi-page PDF ───────────────────────────
message("Generating plots (LOESS smoothing may take a minute)...")

plots <- lapply(seq_len(nrow(param_meta)), function(i) {
  message("  Plotting: ", param_meta$label[i])
  plot_param(param_meta[i, ])
})
plots <- Filter(Negate(is.null), plots)
message("Plots built: ", length(plots))

n_per_page <- 4
n_pages    <- ceiling(length(plots) / n_per_page)

out_pdf <- "C:/Users/fc809/Downloads/NHANES_CBC_age_trends.pdf"
pdf(out_pdf, width = 11, height = 9)

for (pg in seq_len(n_pages)) {
  idx   <- ((pg - 1) * n_per_page + 1) : min(pg * n_per_page, length(plots))
  p_sub <- plots[idx]
  while (length(p_sub) < n_per_page) p_sub <- c(p_sub, list(plot_spacer()))
  print(
    wrap_plots(p_sub, ncol = 2) +
      plot_annotation(
        title    = "NHANES CBC Parameters vs Age  |  Weighted LOESS ± IQR",
        subtitle = paste0("Cycles 2005–2023  |  Excluded: cancer, anaemia, renal failure, liver disease, autoimmune, HIV, ",
                          "pregnancy, CRP>10, smoking, BMI≥35, heavy alcohol, diabetes  |  Blue = Male, Red = Female"),
        caption  = paste0("Exclusions: cancer (MCQ220), anaemia (MCQ092), liver disease (MCQ160L), autoimmune (MCQ160A/N), ",
                          "renal failure (KIQ022/025), HIV (LBXHIV), pregnancy (RHD143/RHQ141), CRP>10mg/L, ",
                          "smoking (SMQ020), BMI≥35, heavy alcohol (>14 drinks/wk), diabetes (DIQ010). ",
                          "Tukey 3×IQR outlier fence applied per parameter."),
        theme = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, colour = "grey40"),
          plot.caption  = element_text(size = 7, colour = "grey50")
        )
      )
  )
}

dev.off()
message("Done! PDF saved to: ", out_pdf)