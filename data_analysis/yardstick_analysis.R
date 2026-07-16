library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(readxl)
library(readr)
library(tidyverse)
library(bayesplot)
library(future)
library(furrr)
library(car)
library(lme4)        # mixed-effects sensitivity model
library(lmtest)      # Breusch-Pagan heteroscedasticity test
library(segmented)   # piecewise regression for time-dependence test
library(DHARMa)      # simulation-based residual + formal outlier tests

# ── CONFIGURATION ─────────────────────────────────────────────────────────────

MIN_QUERY_ALIGN     <- 80          # minimum query alignment %
MIN_REF_ALIGN       <- 80          # minimum reference alignment %
MIN_FAMILY_PAIRS    <- 10          # minimum pairs to attempt a family yardstick
MIN_CLASS_PAIRS     <- 20          # minimum pairs to attempt a class yardstick
GENUS_CAP           <- 25          # max pairs per genus — matches PMM §2.5.2
COMMON_WINDOW       <- c(5, 15)    # MYA — overlapping support across classes
R2_THRESHOLD_FAMILY <- 0.30        # families below this get LMM intercept only
# How the family "time-informative" cutoff is chosen (Sundaram review):
#   "fixed"    — use R2_THRESHOLD_FAMILY as-is (0.30). The quantile of 0.30
#                within the observed R² distribution is still reported for
#                justification.
#   "quantile" — set the cutoff to the R2_QUANTILE_TARGET quantile of the
#                observed family R² distribution (data-driven justification).
R2_CUTOFF_MODE      <- "fixed"     # "fixed" or "quantile"
R2_QUANTILE_TARGET  <- 0.50        # used when R2_CUTOFF_MODE == "quantile"
PI_LEVELS           <- c(0.95, 0.99)  # PI levels for outlier flagging
FDR_ALPHA           <- 0.05        # Benjamini–Hochberg threshold
REF_TIMES_MYA       <- c(1, 2, 5, 10, 15)  # reference points for class comparison
SEED                <- 12345

OUT_DIR <- "yardstick_v4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "diagnostics"), showWarnings = FALSE)

# ── DATA LOAD AND CLEANING ────────────────────────────────────────────────────

setwd("/scratch/gautschi/allen715/GD_models")
traits     <- read_excel("/scratch/gautschi/allen715/GD_models/final_dataset/Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
divergence <- read_csv("/scratch/gautschi/allen715/GD_models/final_dataset/Complete_divergence_05-14-26.csv")

colnames(traits)     <- make.names(colnames(traits))
colnames(divergence) <- make.names(colnames(divergence))

cat(sprintf("Loaded %d species with traits\n",        nrow(traits)))
cat(sprintf("Loaded %d species pairs with divergence data\n", nrow(divergence)))

# 1) parse genome metrics --------------------------------------------------

parse_n50_bp <- function(x) {
  if (is.na(x)) return(NA_real_)
  x     <- str_trim(x)
  value <- as.numeric(str_extract(x, "[0-9.]+"))
  unit  <- str_extract(tolower(x), "kb|mb|gb")
  if (is.na(value) | is.na(unit)) return(NA_real_)
  value * case_when(
    unit == "kb" ~ 1e3,
    unit == "mb" ~ 1e6,
    unit == "gb" ~ 1e9,
    TRUE         ~ 1
  )
}

traits <- traits %>%
  mutate(
    Genome.size           = suppressWarnings(as.numeric(Genome.size)),
    Genome.repeat.content = suppressWarnings(as.numeric(Genome.repeat.content)),
    Genome.repeat.content = na_if(Genome.repeat.content, 0),
    Contig_N50_bp         = map_dbl(Genome.Contig.N50,   parse_n50_bp),
    Scaffold_N50_bp       = map_dbl(Genome.Scaffold.N50, parse_n50_bp),
    N50_bp                = coalesce(Scaffold_N50_bp, Contig_N50_bp)
  )

# 2) exclude flagged species -----------------------------------------------

excluded_species <- traits %>%
  filter(Exclude. %in% c("yes", "Yes")) %>%
  pull(accession)

# Note: prior outputs flagged Bos gaurus / B. frontalis as a yardstick outlier;
# B. frontalis is the domestic form of B. gaurus and the TimeTree estimate of
# 0.21 MYA reflects domestication, not species divergence. If not already in
# Exclude. = yes, this pair will reappear as an outlier and should be reviewed.
flag_domestic_pair <- function(div) {
  hits <- div %>%
    filter(
      (str_detect(species1, "frontalis") & str_detect(species2, "gaurus")) |
      (str_detect(species1, "gaurus") & str_detect(species2, "frontalis"))
    )
  if (nrow(hits) > 0) {
    cat("  WARNING: B. gaurus / B. frontalis pair is in the dataset.\n")
    cat("    frontalis is the domestic form of gaurus and should likely be\n")
    cat("    flagged in the Exclude. column of the traits file.\n")
  }
}
flag_domestic_pair(divergence)

cat(sprintf("Excluding %d flagged species\n", length(excluded_species)))

divergence_clean <- divergence %>%
  filter(!genome1 %in% excluded_species,
         !genome2 %in% excluded_species)

cat(sprintf("Pairs after species exclusion: %d\n", nrow(divergence_clean)))

# 3) alignment quality filter ----------------------------------------------

divergence_clean <- divergence_clean %>%
  filter(
    !is.na(query_alignment_percent),
    !is.na(ref_alignment_percent),
    query_alignment_percent >= MIN_QUERY_ALIGN,
    ref_alignment_percent   >= MIN_REF_ALIGN
  )

cat(sprintf("Pairs after alignment filters (>=%d%% both genomes): %d\n",
            MIN_QUERY_ALIGN, nrow(divergence_clean)))

# 4) name normalisation ----------------------------------------------------

normalize_species <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("[^A-Za-z ]", "") %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

traits <- traits %>%
  mutate(species_normalized = normalize_species(species),
         genus_normalized   = str_to_lower(genus))

divergence_clean <- divergence_clean %>%
  mutate(
    species1_norm = normalize_species(species1),
    species2_norm = normalize_species(species2)
  )

# 5) randomize sp1 / sp2 assignment ----------------------------------------

set.seed(SEED)

divergence_clean <- divergence_clean %>%
  rowwise() %>%
  mutate(
    swap       = sample(c(TRUE, FALSE), 1),
    sp1        = ifelse(swap, species2_norm, species1_norm),
    sp2        = ifelse(swap, species1_norm, species2_norm),
    genome_sp1 = ifelse(swap, genome2, genome1),
    genome_sp2 = ifelse(swap, genome1, genome2)
  ) %>%
  ungroup() %>%
  dplyr::select(-swap)

# 6) merge traits to pairs -------------------------------------------------

trait_cols <- c(
  "Body.size..kg.", "Generation.time..years.", "Clutch.size",
  "Population.Density", "Genome.size", "Genome.repeat.content",
  "Genome.GC.content", "area_km2", "N50_bp"
)

traits_sp1 <- traits %>%
  dplyr::select(accession, Class, Order, Family, genus, all_of(trait_cols)) %>%
  rename(Class_sp1 = Class, Order_sp1 = Order, Family_sp1 = Family, genus_sp1 = genus) %>%
  rename_with(~ paste0(.x, "_sp1"), all_of(trait_cols))

traits_sp2 <- traits %>%
  dplyr::select(accession, Class, Order, Family, genus, all_of(trait_cols)) %>%
  rename(Class_sp2 = Class, Order_sp2 = Order, Family_sp2 = Family, genus_sp2 = genus) %>%
  rename_with(~ paste0(.x, "_sp2"), all_of(trait_cols))

merged_data <- divergence_clean %>%
  left_join(traits_sp1, by = c("genome_sp1" = "accession")) %>%
  left_join(traits_sp2, by = c("genome_sp2" = "accession")) %>%
  mutate(
    Class_final  = coalesce(Class_sp1,  Class_sp2),
    Order_final  = coalesce(Order_sp1,  Order_sp2),
    Family_final = coalesce(Family_sp1, Family_sp2),
    genus_final  = str_to_lower(coalesce(genus_sp1, genus_sp2))
  )

# pair-level trait averages and differences (used by downstream code only)
for (col in trait_cols) {
  sp1_col <- paste0(col, "_sp1")
  sp2_col <- paste0(col, "_sp2")
  v1 <- as.numeric(merged_data[[sp1_col]])
  v2 <- as.numeric(merged_data[[sp2_col]])
  merged_data[[paste0(col, "_avg")]]  <- ifelse(
    is.na(v1) & is.na(v2), NA, rowMeans(cbind(v1, v2), na.rm = TRUE))
  merged_data[[paste0(col, "_diff")]] <- ifelse(
    is.na(v1) | is.na(v2), NA, abs(v1 - v2))
}

# ── BUILD YARDSTICK DATASET ───────────────────────────────────────────────────
# All four classes, no divergence-time cap. Alignment filters retained.
# k2p > 0 and timetree_div > 0 enforced because log(0) = -Inf silently empties
# the regression (this was the reason Amphibia disappeared from earlier outputs).

yardstick_data <- merged_data %>%
  filter(
    Class_final %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"),
    !is.na(k2p),
    !is.na(timetree_div),
    k2p > 0,
    timetree_div > 0
  ) %>%
  mutate(
    log_k2p      = log(k2p),
    log_div_time = log(timetree_div),
    Class_final  = factor(Class_final,
                          levels = c("Aves", "Mammalia", "Reptilia", "Amphibia"))
  )

cat("\n=== YARDSTICK DATASET (FULL) ===\n")
yardstick_data %>%
  group_by(Class_final) %>%
  summarize(
    n_pairs        = n(),
    n_genera       = n_distinct(genus_final),
    median_divtime = median(timetree_div),
    range_divtime  = sprintf("%.2f – %.2f My",
                             min(timetree_div), max(timetree_div)),
    .groups = "drop"
  ) %>%
  print()

# ── APPLY GENUS CAP (PRIMARY DATASET) ─────────────────────────────────────────
# Matches the PMM cap (§2.5.2). Falconidae (Falco) and Psittacidae (Pyrrhura)
# previously dominated Aves; the cap brings any genus to ≤25 pairs.

set.seed(SEED)
yardstick_data_capped <- yardstick_data %>%
  group_by(Class_final, genus_final) %>%
  slice_sample(n = GENUS_CAP, replace = FALSE) %>%
  ungroup()

# log what the cap did
genus_cap_summary <- yardstick_data %>%
  group_by(Class_final, genus_final) %>%
  summarize(n_before = n(), .groups = "drop") %>%
  mutate(n_after = pmin(n_before, GENUS_CAP),
         n_dropped = n_before - n_after) %>%
  filter(n_dropped > 0) %>%
  arrange(desc(n_dropped))

cat(sprintf("\n=== GENUS CAP (max %d pairs/genus) ===\n", GENUS_CAP))
cat(sprintf("Genera capped: %d\n", nrow(genus_cap_summary)))
if (nrow(genus_cap_summary) > 0) {
  cat("Top genera capped (n_before -> n_after):\n")
  print(genus_cap_summary %>% head(10), n = 10)
}

cat("\n=== YARDSTICK DATASET (GENUS-CAPPED, PRIMARY) ===\n")
yardstick_data_capped %>%
  group_by(Class_final) %>%
  summarize(
    n_pairs  = n(),
    n_genera = n_distinct(genus_final),
    range_divtime = sprintf("%.2f – %.2f My",
                            min(timetree_div), max(timetree_div)),
    .groups = "drop"
  ) %>%
  print()

write_csv(genus_cap_summary, file.path(OUT_DIR, "diagnostics", "genus_cap_summary.csv"))

# ── STEP 1: FIT CLASS-LEVEL REGRESSIONS ───────────────────────────────────────
# OLS log(k2p) ~ log(divergence_time) per class, on the GENUS-CAPPED dataset.
# Sensitivity fits (uncapped, common-window, LMM) follow in Step 1b.

fit_class_regression <- function(df, class_name, label = "") {
  df_class <- df %>% filter(Class_final == class_name)

  if (nrow(df_class) < MIN_CLASS_PAIRS) {
    cat(sprintf("  %s [%s]: insufficient pairs (%d) — skipping\n",
                class_name, label, nrow(df_class)))
    return(NULL)
  }

  fit <- lm(log_k2p ~ log_div_time, data = df_class)

  cat(sprintf("\n=== CLASS REGRESSION: %s [%s] (n = %d) ===\n",
              class_name, label, nrow(df_class)))
  cat(sprintf("  Intercept:    %.3f\n", coef(fit)[1]))
  cat(sprintf("  log_div_time: %.3f\n", coef(fit)[2]))
  cat(sprintf("  R-squared:    %.3f\n", summary(fit)$r.squared))
  cat(sprintf("  Residual SD:  %.3f\n", summary(fit)$sigma))

  if (nrow(df_class) < 100) {
    cat(sprintf("  NOTE: n = %d is marginal; interpret PIs with caution.\n",
                nrow(df_class)))
  }

  return(fit)
}

cat("\n=== PRIMARY (GENUS-CAPPED) CLASS REGRESSIONS ===")

class_fits <- map(
  levels(yardstick_data_capped$Class_final),
  ~ fit_class_regression(yardstick_data_capped, .x, label = "primary, capped")
) %>%
  setNames(levels(yardstick_data_capped$Class_final))

# per-class divergence time ranges (for extrapolation warnings in predictions)
class_div_ranges <- yardstick_data_capped %>%
  group_by(Class_final) %>%
  summarize(
    div_min = min(timetree_div, na.rm = TRUE),
    div_max = max(timetree_div, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  { setNames(
      map2(.$div_min, .$div_max, ~ list(min = .x, max = .y)),
      .$Class_final
    ) }

# ── STEP 1b: SENSITIVITY FITS ─────────────────────────────────────────────────

# Sensitivity 1: uncapped OLS (prior method)
cat("\n=== SENSITIVITY 1: UNCAPPED OLS ===")
class_fits_uncapped <- map(
  levels(yardstick_data$Class_final),
  ~ fit_class_regression(yardstick_data, .x, label = "uncapped")
) %>%
  setNames(levels(yardstick_data$Class_final))

# Sensitivity 2: common divergence-time window (5–15 MYA)
# The only honest way to compare slopes across classes, given disparate ranges.
cat(sprintf("\n=== SENSITIVITY 2: COMMON WINDOW (%g–%g MYA) OLS ===",
            COMMON_WINDOW[1], COMMON_WINDOW[2]))
common_window_data <- yardstick_data %>%
  filter(timetree_div >= COMMON_WINDOW[1],
         timetree_div <= COMMON_WINDOW[2])

cat(sprintf("\nPairs in %g–%g MYA window by class:\n",
            COMMON_WINDOW[1], COMMON_WINDOW[2]))
print(table(common_window_data$Class_final))

class_fits_common <- map(
  levels(common_window_data$Class_final),
  ~ fit_class_regression(common_window_data, .x, label = "common window")
) %>%
  setNames(levels(common_window_data$Class_final))

# Sensitivity 3: LMM with genus random intercept (full data, all pairs)
# Compared to OLS, this absorbs genus-level rate variation into a random
# intercept so the fixed-effect slope reflects only within-genus structure
# inflated by between-genus structure to the extent it varies systematically
# with divergence time.
cat("\n=== SENSITIVITY 3: LMM (RANDOM INTERCEPT BY GENUS) ===")
fit_class_lmm <- function(df, class_name) {
  df_class <- df %>% filter(Class_final == class_name)
  if (nrow(df_class) < MIN_CLASS_PAIRS) return(NULL)
  # require at least 5 genera for a meaningful random effect
  if (n_distinct(df_class$genus_final) < 5) {
    cat(sprintf("  %s: too few genera (%d) for LMM\n",
                class_name, n_distinct(df_class$genus_final)))
    return(NULL)
  }
  fit <- tryCatch(
    lmer(log_k2p ~ log_div_time + (1 | genus_final), data = df_class,
         REML = TRUE),
    error = function(e) { cat("  LMM error: ", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)

  fe <- fixef(fit)
  cat(sprintf("\n%s LMM (n = %d, n_genera = %d):\n",
              class_name, nrow(df_class), n_distinct(df_class$genus_final)))
  cat(sprintf("  Fixed intercept: %.3f\n", fe[1]))
  cat(sprintf("  Fixed slope:     %.3f\n", fe[2]))
  cat(sprintf("  Genus SD:        %.3f\n",
              attr(VarCorr(fit)$genus_final, "stddev")))
  cat(sprintf("  Residual SD:     %.3f\n", sigma(fit)))
  return(fit)
}

class_fits_lmm <- map(
  levels(yardstick_data$Class_final),
  ~ fit_class_lmm(yardstick_data, .x)
) %>%
  setNames(levels(yardstick_data$Class_final))

# Pooled all-class fit (for Amphibia and pan-vertebrate baseline)
cat("\n=== POOLED (ALL-CLASS) REGRESSION ===\n")
pooled_fit   <- lm(log_k2p ~ log_div_time, data = yardstick_data_capped)
pooled_sigma <- summary(pooled_fit)$sigma
cat(sprintf("  n pairs:    %d (capped)\n", nrow(yardstick_data_capped)))
cat(sprintf("  Intercept:  %.3f\n", coef(pooled_fit)[1]))
cat(sprintf("  Slope:      %.3f\n", coef(pooled_fit)[2]))
cat(sprintf("  R-squared:  %.3f\n", summary(pooled_fit)$r.squared))
cat(sprintf("  Residual SD: %.3f\n", pooled_sigma))

pooled_div_range <- list(
  min = min(yardstick_data_capped$timetree_div, na.rm = TRUE),
  max = max(yardstick_data_capped$timetree_div, na.rm = TRUE)
)

# ── STEP 1c: SLOPE / INTERCEPT COMPARISON TABLE ───────────────────────────────
# Side-by-side comparison of all four fits, per class. The Discussion can
# cite which estimates are stable across methods (most credible) vs. which
# are sensitive to method choice (interpret cautiously).

extract_lm_params <- function(fit, method, class_name) {
  if (is.null(fit)) {
    return(tibble(class = class_name, method = method,
                  n = NA_integer_, intercept = NA_real_, slope = NA_real_,
                  slope_se = NA_real_, r_squared = NA_real_, residual_sd = NA_real_))
  }
  s <- summary(fit)
  tibble(
    class       = class_name,
    method      = method,
    n           = length(residuals(fit)),
    intercept   = unname(coef(fit)[1]),
    slope       = unname(coef(fit)[2]),
    slope_se    = unname(coef(s)[2, "Std. Error"]),
    r_squared   = s$r.squared,
    residual_sd = s$sigma
  )
}

extract_lmm_params <- function(fit, method, class_name) {
  if (is.null(fit)) {
    return(tibble(class = class_name, method = method,
                  n = NA_integer_, intercept = NA_real_, slope = NA_real_,
                  slope_se = NA_real_, r_squared = NA_real_, residual_sd = NA_real_))
  }
  fe <- fixef(fit)
  se <- summary(fit)$coefficients[2, "Std. Error"]
  tibble(
    class       = class_name,
    method      = method,
    n           = nobs(fit),
    intercept   = unname(fe[1]),
    slope       = unname(fe[2]),
    slope_se    = se,
    r_squared   = NA_real_,    # LMM marginal R² requires MuMIn; omitted for portability
    residual_sd = sigma(fit)
  )
}

slope_comparison <- bind_rows(
  map_dfr(names(class_fits),           ~ extract_lm_params(class_fits[[.x]],          "primary_capped",  .x)),
  map_dfr(names(class_fits_uncapped),  ~ extract_lm_params(class_fits_uncapped[[.x]], "uncapped",        .x)),
  map_dfr(names(class_fits_common),    ~ extract_lm_params(class_fits_common[[.x]],   "common_window",   .x)),
  map_dfr(names(class_fits_lmm),       ~ extract_lmm_params(class_fits_lmm[[.x]],     "lmm_genus_ri",    .x))
) %>%
  arrange(class, method)

cat("\n=== SLOPE / INTERCEPT COMPARISON ACROSS METHODS ===\n")
print(slope_comparison, n = 30)
write_csv(slope_comparison, file.path(OUT_DIR, "slope_comparison_methods.csv"))

# ── STEP 2: RESIDUAL DIAGNOSTICS (PRIMARY FITS) ───────────────────────────────
# Normality, heteroscedasticity, and time-dependence checks per class.

cat("\n=== RESIDUAL DIAGNOSTICS (PRIMARY FITS) ===\n")

diagnose_class <- function(fit, df_class, class_name) {
  if (is.null(fit)) return(NULL)
  resid_vec <- residuals(fit)
  fitted_vec <- fitted(fit)

  # Normality
  sw <- if (length(resid_vec) >= 3 && length(resid_vec) <= 5000) {
    shapiro.test(resid_vec)
  } else NULL

  # Heteroscedasticity (Breusch-Pagan)
  bp <- tryCatch(bptest(fit), error = function(e) NULL)

  # Time-dependence: segmented regression
  seg <- tryCatch(
    segmented(fit, seg.Z = ~ log_div_time),
    error = function(e) NULL
  )
  seg_breakpoint_mya <- if (!is.null(seg) && !is.null(seg$psi)) {
    exp(seg$psi[1, "Est."])
  } else NA_real_
  # Davies test for non-constant slope (more portable than pscore.test
  # across segmented package versions)
  seg_pvalue <- tryCatch(
    davies.test(fit, seg.Z = ~ log_div_time)$p.value,
    error = function(e) NA_real_
  )

  # Residual SD across binned fitted values (visual heteroscedasticity)
  bin_summary <- tibble(fitted = fitted_vec, residual = resid_vec) %>%
    mutate(bin = cut(fitted, breaks = quantile(fitted, probs = seq(0, 1, 0.2),
                                                na.rm = TRUE),
                    include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarize(n = n(), bin_sd = sd(residual), .groups = "drop")

  cat(sprintf("\n%s:\n", class_name))
  if (!is.null(sw))  cat(sprintf("  Shapiro–Wilk W = %.3f, p = %.4f\n",
                                 sw$statistic, sw$p.value))
  if (!is.null(bp))  cat(sprintf("  Breusch–Pagan BP = %.3f, df = %d, p = %.4f\n",
                                 bp$statistic, bp$parameter, bp$p.value))
  if (!is.na(seg_breakpoint_mya)) {
    cat(sprintf("  Segmented regression breakpoint: %.2f MYA (Davies p = %.4f)\n",
                seg_breakpoint_mya, seg_pvalue))
  } else {
    cat("  Segmented regression: no significant breakpoint detected\n")
  }
  cat("  Residual SD by binned fitted value (heteroscedasticity check):\n")
  print(bin_summary)

  list(
    class      = class_name,
    shapiro_W  = if (!is.null(sw)) unname(sw$statistic) else NA_real_,
    shapiro_p  = if (!is.null(sw)) sw$p.value           else NA_real_,
    bp_stat    = if (!is.null(bp)) unname(bp$statistic) else NA_real_,
    bp_p       = if (!is.null(bp)) bp$p.value           else NA_real_,
    seg_breakpoint_mya = seg_breakpoint_mya,
    seg_pvalue = seg_pvalue,
    bin_summary = bin_summary
  )
}

diagnostics <- map2(
  class_fits,
  names(class_fits),
  ~ diagnose_class(.x, yardstick_data_capped %>% filter(Class_final == .y), .y)
)

diagnostics_table <- map_dfr(
  diagnostics[!sapply(diagnostics, is.null)],
  ~ tibble(class = .x$class,
           shapiro_W = .x$shapiro_W, shapiro_p = .x$shapiro_p,
           bp_stat = .x$bp_stat, bp_p = .x$bp_p,
           seg_breakpoint_mya = .x$seg_breakpoint_mya,
           seg_pvalue = .x$seg_pvalue)
)
write_csv(diagnostics_table, file.path(OUT_DIR, "diagnostics", "residual_diagnostics.csv"))

# Q-Q plots saved as supplementary diagnostic figures
qq_plots <- imap(class_fits, function(fit, class_name) {
  if (is.null(fit)) return(NULL)
  df <- tibble(resid = residuals(fit))
  ggplot(df, aes(sample = resid)) +
    stat_qq(alpha = 0.5) + stat_qq_line(color = "red") +
    labs(title = sprintf("Q–Q plot of residuals: %s", class_name),
         x = "Theoretical quantile", y = "Sample residual") +
    theme_minimal(base_size = 12)
})

walk2(qq_plots, names(qq_plots), function(p, cn) {
  if (!is.null(p)) {
    ggsave(file.path(OUT_DIR, "diagnostics", sprintf("qq_%s.png", cn)),
           p, width = 5, height = 5, dpi = 150)
  }
})

# ── STEP 3: COMPUTE YARDSTICK (PRIMARY FITS, FULL DATA) ───────────────────────
# Even though fits are from genus-capped data, residuals and outlier flags
# are applied to ALL pairs in yardstick_data — the cap controls the fit, not
# the inferences about individual pairs.

compute_yardstick <- function(df, fit, class_name) {
  df_class <- df %>% filter(Class_final == class_name)
  if (is.null(fit) || nrow(df_class) == 0) return(NULL)

  sigma <- summary(fit)$sigma
  df_class <- df_class %>%
    mutate(
      fitted_log_k2p = predict(fit, newdata = df_class),
      residual       = log_k2p - fitted_log_k2p,
      std_residual   = residual / sigma
    )

  resid_summary <- tibble(
    class          = class_name,
    n_pairs_fit    = length(residuals(fit)),
    n_pairs_full   = nrow(df_class),
    resid_mean     = mean(df_class$residual),
    resid_sd       = sd(df_class$residual),
    resid_p025     = quantile(df_class$residual, 0.025),
    resid_p25      = quantile(df_class$residual, 0.25),
    resid_p75      = quantile(df_class$residual, 0.75),
    resid_p975     = quantile(df_class$residual, 0.975),
    std_resid_p025 = quantile(df_class$std_residual, 0.025),
    std_resid_p975 = quantile(df_class$std_residual, 0.975),
    intercept      = coef(fit)[1],
    slope          = coef(fit)[2],
    r_squared      = summary(fit)$r.squared,
    residual_sd    = sigma,
    div_time_min   = min(df_class$timetree_div),
    div_time_max   = max(df_class$timetree_div)
  )

  list(data = df_class, summary = resid_summary)
}

yardstick_results <- map2(
  class_fits,
  names(class_fits),
  ~ compute_yardstick(yardstick_data, .x, .y)
)

yardstick_table <- bind_rows(map(yardstick_results, ~ .x$summary))
cat("\n=== CLASS-LEVEL YARDSTICK TABLE (PRIMARY FITS) ===\n")
print(yardstick_table)
write_csv(yardstick_table, file.path(OUT_DIR, "class_yardstick_table.csv"))

# ── STEP 4: PREDICT YARDSTICK FOR A NEW PAIR ──────────────────────────────────
# Uses predict(..., interval = "prediction") which computes a proper
# x-dependent PI — wider at the extremes of the divergence-time range,
# narrower near the mean — addressing the prior flat-band concern.

predict_yardstick <- function(div_time_mya, class_name,
                              class_fits, class_div_ranges,
                              alpha = 0.05) {
  fit <- class_fits[[class_name]]
  if (is.null(fit)) {
    cat(sprintf("  No model available for %s\n", class_name))
    return(NULL)
  }

  rng <- class_div_ranges[[class_name]]
  extrap <- !is.null(rng) && (div_time_mya < rng$min || div_time_mya > rng$max)
  if (extrap) {
    cat(sprintf(
      "  WARNING: %.2f My is outside the %s training range (%.2f – %.2f My).\n",
      div_time_mya, class_name, rng$min, rng$max))
  }

  newdata <- tibble(log_div_time = log(div_time_mya))
  pred    <- predict(fit, newdata = newdata,
                     interval = "prediction", level = 1 - alpha)
  tibble(
    class        = class_name,
    div_time_mya = div_time_mya,
    expected_k2p = exp(pred[, "fit"]),
    lower        = exp(pred[, "lwr"]),
    upper        = exp(pred[, "upr"]),
    expected_log = pred[, "fit"],
    lower_log    = pred[, "lwr"],
    upper_log    = pred[, "upr"],
    extrapolation = extrap
  )
}

# ── STEP 5: CROSS-CLASS REFERENCE PREDICTIONS ─────────────────────────────────
# Predicted k2p at standardized divergence times for direct cross-class
# comparison — addresses Goal #2 (compare/contrast with Johns & Avise 1998).
# At each reference time we report the central estimate and the PI; the
# extrapolation flag warns when the reference time is outside a given
# class's training range.

cat("\n=== CROSS-CLASS REFERENCE PREDICTIONS (for Johns & Avise comparison) ===\n")

reference_predictions <- expand_grid(
  div_time_mya = REF_TIMES_MYA,
  class_name   = names(class_fits)
) %>%
  pmap_dfr(function(div_time_mya, class_name) {
    if (is.null(class_fits[[class_name]])) return(NULL)
    predict_yardstick(div_time_mya, class_name, class_fits, class_div_ranges)
  }) %>%
  mutate(
    k2p_pretty = sprintf("%.4f [%.4f – %.4f]", expected_k2p, lower, upper)
  )

# Wide format: rows = reference times, cols = classes
ref_wide <- reference_predictions %>%
  dplyr::select(div_time_mya, class, k2p_pretty) %>%
  pivot_wider(names_from = class, values_from = k2p_pretty)

cat("Predicted k2p [95% PI] by class at standard divergence times:\n")
print(ref_wide)

write_csv(reference_predictions, file.path(OUT_DIR, "reference_predictions_by_class.csv"))
write_csv(ref_wide, file.path(OUT_DIR, "reference_predictions_by_class_wide.csv"))

# ── STEP 6: OUTLIER FLAGGING — 95%, 99%, AND FDR-ADJUSTED ─────────────────────
# Three tiers to distinguish chance-expected flags from substantive findings.

flag_outliers_extended <- function(df, class_fits) {
  df %>%
    group_by(Class_final) %>%
    group_modify(~ {
      cn  <- as.character(unique(.y$Class_final))
      fit <- class_fits[[cn]]
      if (is.null(fit)) return(.x)

      sigma  <- summary(fit)$sigma
      df_fit <- length(residuals(fit)) - 2  # for two-sided t-test

      pred95 <- predict(fit, newdata = .x, interval = "prediction", level = 0.95)
      pred99 <- predict(fit, newdata = .x, interval = "prediction", level = 0.99)

      residual     <- .x$log_k2p - pred95[, "fit"]
      std_residual <- residual / sigma
      p_two_sided  <- 2 * pt(abs(std_residual), df = df_fit, lower.tail = FALSE)
      q_fdr        <- p.adjust(p_two_sided, method = "BH")

      .x %>% mutate(
        expected_log_k2p = pred95[, "fit"],
        pi_lower_95      = pred95[, "lwr"],
        pi_upper_95      = pred95[, "upr"],
        pi_lower_99      = pred99[, "lwr"],
        pi_upper_99      = pred99[, "upr"],
        residual         = residual,
        std_residual     = std_residual,
        p_value          = p_two_sided,
        q_fdr            = q_fdr,
        # Standardized-residual outliers (primary reported method, Sundaram
        # review): |z| > 1.96 (~95%) and |z| > 2.576 (~99%). These use the
        # class residual SD directly and are cross-class comparable.
        outlier_sr_95_low  = std_residual < -1.96,
        outlier_sr_95_high = std_residual >  1.96,
        outlier_sr_95      = abs(std_residual) > 1.96,
        outlier_sr_99_low  = std_residual < -2.576,
        outlier_sr_99_high = std_residual >  2.576,
        outlier_sr_99      = abs(std_residual) > 2.576,
        outlier_95_low   = log_k2p < pi_lower_95,
        outlier_95_high  = log_k2p > pi_upper_95,
        outlier_95       = outlier_95_low | outlier_95_high,
        outlier_99_low   = log_k2p < pi_lower_99,
        outlier_99_high  = log_k2p > pi_upper_99,
        outlier_99       = outlier_99_low | outlier_99_high,
        outlier_fdr      = q_fdr < FDR_ALPHA
      )
    }) %>%
    ungroup()
}

yardstick_flagged <- flag_outliers_extended(yardstick_data, class_fits)

# Pooled-model flags for Amphibia (no class-specific model)
pooled_pred <- predict(pooled_fit, newdata = yardstick_flagged,
                       interval = "prediction", level = 0.95)
yardstick_flagged <- yardstick_flagged %>%
  mutate(
    pooled_expected_log = pooled_pred[, "fit"],
    pooled_pi_lower     = pooled_pred[, "lwr"],
    pooled_pi_upper     = pooled_pred[, "upr"],
    pooled_residual     = log_k2p - pooled_pred[, "fit"],
    pooled_std_residual = pooled_residual / pooled_sigma,
    pooled_outlier_low  = log_k2p < pooled_pred[, "lwr"],
    pooled_outlier_high = log_k2p > pooled_pred[, "upr"],
    pooled_outlier      = pooled_outlier_low | pooled_outlier_high
  )

outlier_summary <- yardstick_flagged %>%
  group_by(Class_final) %>%
  summarize(
    n_pairs           = n(),
    n_outlier_sr_95   = sum(outlier_sr_95, na.rm = TRUE),
    n_outlier_sr_99   = sum(outlier_sr_99, na.rm = TRUE),
    n_outlier_95      = sum(outlier_95,  na.rm = TRUE),
    n_outlier_99      = sum(outlier_99,  na.rm = TRUE),
    n_outlier_fdr     = sum(outlier_fdr, na.rm = TRUE),
    pct_outlier_sr_95 = 100 * mean(outlier_sr_95, na.rm = TRUE),
    pct_outlier_sr_99 = 100 * mean(outlier_sr_99, na.rm = TRUE),
    pct_outlier_95    = 100 * mean(outlier_95,  na.rm = TRUE),
    pct_outlier_99    = 100 * mean(outlier_99,  na.rm = TRUE),
    pct_outlier_fdr   = 100 * mean(outlier_fdr, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== OUTLIER SUMMARY (STANDARDIZED RESIDUAL + PI + FDR) ===\n")
cat("SR |z|>1.96 (primary): expected ~5% by chance | SR |z|>2.576: expected ~1%\n")
cat("95% PI: expected ~5% by chance | 99% PI: expected ~1% by chance | FDR: BH at q<0.05\n")
print(outlier_summary)
write_csv(outlier_summary, file.path(OUT_DIR, "outlier_summary_by_class.csv"))

# Flagged pairs CSV — 95% threshold (matches prior output for compatibility)
write_csv(
  yardstick_flagged %>%
    filter(outlier_95) %>%
    dplyr::select(
      Class_final, Family_final, genus_final, sp1, sp2,
      timetree_div, k2p,
      expected_log_k2p, pi_lower_95, pi_upper_95,
      pi_lower_99, pi_upper_99,
      residual, std_residual, p_value, q_fdr,
      outlier_95_low, outlier_95_high,
      outlier_99, outlier_fdr
    ),
  file.path(OUT_DIR, "flagged_outlier_pairs.csv")
)

# 99%-threshold pairs (the substantive set, much smaller)
write_csv(
  yardstick_flagged %>%
    filter(outlier_99) %>%
    dplyr::select(
      Class_final, Family_final, genus_final, sp1, sp2,
      timetree_div, k2p,
      expected_log_k2p, pi_lower_99, pi_upper_99,
      residual, std_residual, p_value, q_fdr,
      outlier_99_low, outlier_99_high
    ),
  file.path(OUT_DIR, "flagged_outlier_pairs_99pct.csv")
)

# FDR-significant pairs
write_csv(
  yardstick_flagged %>%
    filter(outlier_fdr) %>%
    dplyr::select(
      Class_final, Family_final, genus_final, sp1, sp2,
      timetree_div, k2p,
      expected_log_k2p, residual, std_residual, p_value, q_fdr
    ),
  file.path(OUT_DIR, "flagged_outlier_pairs_fdr.csv")
)

# Standardized-residual pairs (primary reported method: |z| > 1.96)
write_csv(
  yardstick_flagged %>%
    filter(outlier_sr_95) %>%
    dplyr::select(
      Class_final, Family_final, genus_final, sp1, sp2,
      timetree_div, k2p,
      expected_log_k2p, residual, std_residual, p_value, q_fdr,
      outlier_sr_95_low, outlier_sr_95_high, outlier_sr_99
    ) %>%
    arrange(desc(abs(std_residual))),
  file.path(OUT_DIR, "flagged_outlier_pairs_stdresid.csv")
)

# ── STEP 6b: DHARMa CROSS-CHECK OF OUTLIERS ───────────────────────────────────
# Sundaram review suggested a simulation-based formal outlier test as an
# independent check on the standardized-residual flags. DHARMa's
# simulateResiduals() rescales residuals to a uniform(0,1) scale via posterior
# predictive simulation, and testOutliers() runs a formal binomial test on the
# number of points falling at the 0/1 boundary. It accepts both lm and lmer
# (merMod) objects, so we run it on (i) the primary genus-capped class OLS fits
# and (ii) the genus-random-intercept LMMs from Step 1b. This does NOT replace
# the standardized-residual flags reported per pair; it provides a class-level
# p-value for "are there more outliers than a correctly-specified model
# predicts?" for the manuscript / a quantitatively-minded reviewer.

cat("\n=== STEP 6b: DHARMa FORMAL OUTLIER TEST (CROSS-CHECK) ===\n")

set.seed(SEED)

run_dharma_outlier_test <- function(fit, class_name, model_label) {
  if (is.null(fit)) return(NULL)
  res <- tryCatch({
    sim <- simulateResiduals(fittedModel = fit, n = 250, plot = FALSE)
    # margin = "both" tests deficits and excesses; binomial test is the default
    ot  <- testOutliers(sim, type = "binomial", plot = FALSE)
    du  <- testUniformity(sim, plot = FALSE)   # KS test on residual uniformity
    tibble(
      class            = class_name,
      model            = model_label,
      n_obs            = length(residuals(fit)),
      outlier_p        = unname(ot$p.value),
      outlier_freq_obs = unname(ot$estimate[1]),
      outlier_test     = ot$method,
      uniformity_ks_p  = unname(du$p.value)
    )
  }, error = function(e) {
    cat(sprintf("  DHARMa failed for %s [%s]: %s\n",
                class_name, model_label, conditionMessage(e)))
    NULL
  })
  if (!is.null(res)) {
    cat(sprintf("  %-10s [%-13s]: testOutliers p = %.4f | uniformity KS p = %.4f\n",
                class_name, model_label, res$outlier_p, res$uniformity_ks_p))
  }
  res
}

dharma_outlier_results <- bind_rows(
  map_dfr(names(class_fits),
          ~ run_dharma_outlier_test(class_fits[[.x]],     .x, "class OLS (capped)")),
  map_dfr(names(class_fits_lmm),
          ~ run_dharma_outlier_test(class_fits_lmm[[.x]], .x, "genus LMM"))
)

if (nrow(dharma_outlier_results) > 0) {
  cat("\nDHARMa outlier-test summary (p < 0.05 => more boundary outliers than expected):\n")
  print(dharma_outlier_results)
  write_csv(dharma_outlier_results,
            file.path(OUT_DIR, "diagnostics", "dharma_outlier_tests.csv"))

  # Concordance note: compare DHARMa's class-level verdict to the SR flag rate.
  sr_rate <- yardstick_flagged %>%
    group_by(Class_final) %>%
    summarize(pct_sr_95 = 100 * mean(outlier_sr_95, na.rm = TRUE), .groups = "drop")
  cat("\nStandardized-residual flag rate by class (compare to DHARMa verdict):\n")
  print(sr_rate)
} else {
  cat("  No DHARMa results produced (all fits NULL?).\n")
}

# ── STEP 7: CLADE-LEVEL OUTLIER SUMMARY ───────────────────────────────────────
# Identify genera with multiple outlier pairs — these are the substantive
# clade-level signals (vs. singleton flags that are likely chance).

genus_outlier_summary <- yardstick_flagged %>%
  group_by(Class_final, Family_final, genus_final) %>%
  summarize(
    n_pairs            = n(),
    n_outliers_95      = sum(outlier_95,      na.rm = TRUE),
    n_outliers_99      = sum(outlier_99,      na.rm = TRUE),
    n_outliers_fdr     = sum(outlier_fdr,     na.rm = TRUE),
    n_outliers_95_low  = sum(outlier_95_low,  na.rm = TRUE),
    n_outliers_95_high = sum(outlier_95_high, na.rm = TRUE),
    mean_std_residual  = mean(std_residual,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_outliers_95 >= 2 | n_outliers_99 >= 1 | n_outliers_fdr >= 1) %>%
  arrange(desc(n_outliers_95), desc(n_outliers_99))

cat("\n=== CLADE-LEVEL OUTLIER SUMMARY (multi-pair outlier genera) ===\n")
print(genus_outlier_summary, n = 50)
write_csv(genus_outlier_summary,
          file.path(OUT_DIR, "genus_outlier_summary.csv"))

# Family-level outlier summary for use as a Discussion anchor
family_outlier_summary <- yardstick_flagged %>%
  group_by(Class_final, Family_final) %>%
  summarize(
    n_pairs            = n(),
    n_outliers_95      = sum(outlier_95,      na.rm = TRUE),
    n_outliers_99      = sum(outlier_99,      na.rm = TRUE),
    n_outliers_95_low  = sum(outlier_95_low,  na.rm = TRUE),
    n_outliers_95_high = sum(outlier_95_high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_outliers_95 >= 2) %>%
  arrange(desc(n_outliers_95))

write_csv(family_outlier_summary,
          file.path(OUT_DIR, "family_outlier_summary.csv"))

# ── STEP 8: FAMILY-LEVEL YARDSTICK (TWO-TIER) ─────────────────────────────────
# Tier A (time-informative): r² ≥ R2_THRESHOLD_FAMILY — full OLS family fit.
# Tier B (intercept-only):  r² < R2_THRESHOLD_FAMILY — fall back to the class
#   slope with a family-specific intercept from a class LMM with family random
#   intercept. This avoids reporting uninterpretable negative-slope family
#   yardsticks (which the prior version produced for Aotidae, Strigidae,
#   Accipitridae, Cervidae, etc.).

fit_family_regression <- function(df, class_name, family_name) {
  df_fam <- df %>%
    filter(Class_final == class_name, Family_final == family_name)
  if (nrow(df_fam) < MIN_FAMILY_PAIRS) return(NULL)
  if (n_distinct(df_fam$timetree_div) < 3) return(NULL)

  fit <- tryCatch(
    lm(log_k2p ~ log_div_time, data = df_fam),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  s <- summary(fit)

  tibble(
    class         = class_name,
    family        = family_name,
    n_pairs       = nrow(df_fam),
    n_genera      = n_distinct(df_fam$genus_final),
    intercept     = coef(fit)[1],
    slope         = coef(fit)[2],
    slope_se      = s$coefficients[2, "Std. Error"],
    r_squared     = s$r.squared,
    residual_sd   = s$sigma,
    resid_p025    = quantile(residuals(fit), 0.025),
    resid_p975    = quantile(residuals(fit), 0.975),
    std_resid_p025 = quantile(residuals(fit) / s$sigma, 0.025),
    std_resid_p975 = quantile(residuals(fit) / s$sigma, 0.975),
    div_time_min  = min(df_fam$timetree_div),
    div_time_max  = max(df_fam$timetree_div)
  )
}

cat("\n=== FAMILY-LEVEL YARDSTICK FITS ===\n")
family_yardsticks_raw <- yardstick_data %>%
  distinct(Class_final, Family_final) %>%
  filter(!is.na(Family_final)) %>%
  pmap_dfr(function(Class_final, Family_final) {
    fit_family_regression(yardstick_data, Class_final, Family_final)
  })

# ── R² CUTOFF JUSTIFICATION (Sundaram review) ─────────────────────────────────
# Rather than assert R² >= 0.30 as an ad hoc "time-informative" threshold, we
# characterize the empirical distribution of family-level R² values and (a)
# report where 0.30 falls as a quantile, and (b) optionally set the cutoff to a
# chosen quantile of that distribution. Either way the choice is now anchored to
# the data instead of a bare number.
r2_vec <- family_yardsticks_raw$r_squared
r2_vec <- r2_vec[is.finite(r2_vec)]

r2_quantile_table <- tibble(
  quantile = c(0.10, 0.25, 0.50, 0.70, 0.75, 0.90),
  r_squared = quantile(r2_vec, probs = c(0.10, 0.25, 0.50, 0.70, 0.75, 0.90),
                       na.rm = TRUE)
)

# Percentile of the fixed 0.30 threshold within the observed distribution:
r2_ecdf              <- ecdf(r2_vec)
r2_threshold_pctile  <- 100 * r2_ecdf(R2_THRESHOLD_FAMILY)
r2_quantile_cutoff   <- as.numeric(quantile(r2_vec, probs = R2_QUANTILE_TARGET,
                                            na.rm = TRUE))

cat("\n=== FAMILY R² DISTRIBUTION (cutoff justification) ===\n")
cat(sprintf("Families with an R² value: %d\n", length(r2_vec)))
cat(sprintf("R² summary: min %.3f | median %.3f | mean %.3f | max %.3f\n",
            min(r2_vec), median(r2_vec), mean(r2_vec), max(r2_vec)))
cat("Quantiles of the family R² distribution:\n")
print(r2_quantile_table)
cat(sprintf("\nThe fixed threshold R² = %.2f sits at the %.1fth percentile of\n",
            R2_THRESHOLD_FAMILY, r2_threshold_pctile))
cat(sprintf("  the observed distribution (i.e. ~%.0f%% of families fall below it).\n",
            r2_threshold_pctile))
cat(sprintf("The %.0fth-percentile (quantile %.2f) cutoff would be R² = %.3f.\n",
            100 * R2_QUANTILE_TARGET, R2_QUANTILE_TARGET, r2_quantile_cutoff))

write_csv(r2_quantile_table,
          file.path(OUT_DIR, "diagnostics", "family_r2_quantiles.csv"))

# Resolve the cutoff actually used, based on R2_CUTOFF_MODE.
r2_cutoff_used <- if (identical(R2_CUTOFF_MODE, "quantile")) {
  r2_quantile_cutoff
} else {
  R2_THRESHOLD_FAMILY
}
cat(sprintf("\nCutoff mode: '%s' -> using R² cutoff = %.3f for tier assignment.\n",
            R2_CUTOFF_MODE, r2_cutoff_used))

family_yardsticks_ols <- family_yardsticks_raw %>%
  mutate(tier = if_else(r_squared >= r2_cutoff_used,
                       "time_informative", "intercept_only"))

cat(sprintf("Families with sufficient data (n >= %d): %d\n",
            MIN_FAMILY_PAIRS, nrow(family_yardsticks_ols)))
cat(sprintf("  Time-informative (r² >= %.2f): %d\n",
            r2_cutoff_used,
            sum(family_yardsticks_ols$tier == "time_informative")))
cat(sprintf("  Intercept-only  (r² < %.2f): %d\n",
            r2_cutoff_used,
            sum(family_yardsticks_ols$tier == "intercept_only")))

# Class LMM with family random intercept — supplies intercept offsets for
# intercept-only families. Returns a per-family BLUP (best linear unbiased
# predictor) of the intercept deviation from the class mean.
fit_class_family_lmm <- function(df, class_name) {
  df_class <- df %>% filter(Class_final == class_name)
  if (nrow(df_class) < MIN_CLASS_PAIRS) return(NULL)
  if (n_distinct(df_class$Family_final) < 3) return(NULL)
  tryCatch(
    lmer(log_k2p ~ log_div_time + (1 | Family_final),
         data = df_class, REML = TRUE),
    error = function(e) NULL
  )
}

class_family_lmms <- map(
  levels(yardstick_data$Class_final),
  ~ fit_class_family_lmm(yardstick_data, .x)
) %>%
  setNames(levels(yardstick_data$Class_final))

extract_family_intercepts <- function(lmm, class_name) {
  if (is.null(lmm)) return(NULL)
  fe   <- fixef(lmm)
  ints <- ranef(lmm)$Family_final
  tibble(
    class            = class_name,
    family           = rownames(ints),
    class_intercept  = unname(fe[1]),
    class_slope      = unname(fe[2]),
    family_offset    = ints[, "(Intercept)"],
    family_intercept = unname(fe[1]) + ints[, "(Intercept)"],
    lmm_residual_sd  = sigma(lmm)
  )
}

family_lmm_intercepts <- map2_dfr(
  class_family_lmms,
  names(class_family_lmms),
  extract_family_intercepts
)

# Combine: time-informative families use OLS fit; intercept-only families use
# class slope + family random intercept.
family_yardsticks <- family_yardsticks_ols %>%
  left_join(family_lmm_intercepts, by = c("class", "family")) %>%
  mutate(
    final_intercept = if_else(tier == "time_informative",
                              intercept, family_intercept),
    final_slope     = if_else(tier == "time_informative",
                              slope, class_slope),
    final_resid_sd  = if_else(tier == "time_informative",
                              residual_sd, lmm_residual_sd)
  )

cat("\nFamily yardstick table (with tier assignment):\n")
print(family_yardsticks %>%
        arrange(class, desc(n_pairs)) %>%
        dplyr::select(class, family, n_pairs, n_genera, tier,
                      slope, slope_se, r_squared,
                      final_intercept, final_slope, final_resid_sd),
      n = 50)

write_csv(family_yardsticks, file.path(OUT_DIR, "family_yardstick_table.csv"))

# ── STEP 9: FIGURES ───────────────────────────────────────────────────────────

class_colors <- c(
  "Aves"     = "#2980B9",
  "Mammalia" = "#C0392B",
  "Reptilia" = "#27AE60",
  "Amphibia" = "#9e22e6"
)

# Figure 1: Class-level yardstick curves with proper x-dependent PIs
plot_divergence_curves <- function(yardstick_data, class_fits) {
  pred_ribbons <- map_dfr(names(class_fits), function(cn) {
    fit <- class_fits[[cn]]
    if (is.null(fit)) return(NULL)
    df_class   <- yardstick_data %>% filter(Class_final == cn)
    time_range <- range(df_class$timetree_div)
    time_seq   <- exp(seq(log(time_range[1]), log(time_range[2]), length.out = 200))
    newdata    <- tibble(log_div_time = log(time_seq))
    pred       <- predict(fit, newdata = newdata, interval = "prediction")
    tibble(
      Class_final = cn,
      div_time    = time_seq,
      fit         = exp(pred[, "fit"]),
      lower       = exp(pred[, "lwr"]),
      upper       = exp(pred[, "upr"])
    )
  })

  ggplot() +
    geom_ribbon(data = pred_ribbons,
                aes(x = div_time, ymin = lower, ymax = upper,
                    fill = Class_final), alpha = 0.15) +
    geom_line(data = pred_ribbons,
              aes(x = div_time, y = fit, color = Class_final),
              linewidth = 1.2) +
    geom_point(data = yardstick_data,
               aes(x = timetree_div, y = k2p, color = Class_final),
               alpha = 0.3, size = 0.8) +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values  = class_colors) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title    = "Genomic divergence yardstick by vertebrate Class",
      subtitle = "Lines = expected k2p; ribbons = 95% prediction interval; primary fit on genus-capped data",
      x        = "Divergence time (My, log scale)",
      y        = "k2p divergence (log scale)",
      color    = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# Figure 2a: raw residual densities per class
plot_residuals_raw <- function(yardstick_results) {
  resid_data <- map_dfr(yardstick_results, ~ {
    if (is.null(.x)) return(NULL)
    .x$data %>% dplyr::select(Class_final, residual)
  })
  ggplot(resid_data, aes(x = residual, fill = Class_final, color = Class_final)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values  = class_colors) +
    labs(title = "Residual distributions by Class (raw log scale)",
         x = "Residual log(k2p)", y = "Density",
         color = NULL, fill = NULL) +
    facet_wrap(~ Class_final, ncol = 2) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

# Figure 2b: standardized residuals — direct cross-class comparison
plot_residuals_standardized <- function(yardstick_results) {
  resid_data <- map_dfr(yardstick_results, ~ {
    if (is.null(.x)) return(NULL)
    .x$data %>% dplyr::select(Class_final, std_residual)
  })
  x_seq <- seq(-4, 4, length.out = 300)
  normal_ref <- tibble(x = x_seq, y = dnorm(x_seq))

  ggplot(resid_data,
         aes(x = std_residual, fill = Class_final, color = Class_final)) +
    geom_density(alpha = 0.20, linewidth = 0.8) +
    geom_line(data = normal_ref, aes(x = x, y = y),
              color = "gray40", linewidth = 0.7, linetype = "dashed",
              inherit.aes = FALSE) +
    geom_vline(xintercept = c(-1.96, 1.96),
               linetype = "dotted", color = "gray50", linewidth = 0.6) +
    geom_vline(xintercept = c(-2.576, 2.576),
               linetype = "dotted", color = "gray30", linewidth = 0.6) +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values  = class_colors) +
    labs(title = "Standardized residuals by Class — cross-class comparable",
         subtitle = "Dashed = N(0,1) reference | dotted gray = ±1.96 (95%) and ±2.576 (99%)",
         x = "Standardized residual (z-score)", y = "Density",
         color = NULL, fill = NULL) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# Figure 3: slope comparison across methods (visual sensitivity check)
plot_slope_comparison <- function(slope_comparison) {
  slope_comparison %>%
    filter(!is.na(slope)) %>%
    ggplot(aes(x = method, y = slope, color = class)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_pointrange(aes(ymin = slope - 1.96 * slope_se,
                        ymax = slope + 1.96 * slope_se),
                    position = position_dodge(width = 0.6),
                    linewidth = 0.7, size = 0.5) +
    scale_color_manual(values = class_colors) +
    labs(title = "Class slope estimates across methods",
         subtitle = "Dashed line at slope = 1 = neutral D = 2kt expectation",
         x = NULL, y = "Slope (log k2p vs log divergence time)",
         color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# Figure 4: reference-time predictions plot (Johns & Avise comparison)
plot_reference_predictions <- function(ref) {
  ggplot(ref %>% filter(!extrapolation),
         aes(x = factor(div_time_mya), y = expected_k2p, color = class)) +
    geom_point(size = 3, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.2, position = position_dodge(width = 0.4),
                  linewidth = 0.7) +
    scale_y_log10() +
    scale_color_manual(values = class_colors) +
    labs(title = "Predicted genomic divergence at standard times",
         subtitle = "For cross-class comparison and Johns & Avise (1998) reference",
         x = "Divergence time (MYA)",
         y = "Predicted k2p [95% PI, log scale]",
         color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

p_curves   <- plot_divergence_curves(yardstick_data, class_fits)
p_residuals <- plot_residuals_raw(yardstick_results)
p_std_resid <- plot_residuals_standardized(yardstick_results)
p_slopes    <- plot_slope_comparison(slope_comparison)
p_refpred   <- plot_reference_predictions(reference_predictions)

ggsave(file.path(OUT_DIR, "divergence_curves_by_class.pdf"),
       p_curves, width = 9, height = 6)
ggsave(file.path(OUT_DIR, "divergence_curves_by_class.png"),
       p_curves, width = 9, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "residual_distributions_by_class.pdf"),
       p_residuals, width = 9, height = 7)
ggsave(file.path(OUT_DIR, "residual_distributions_by_class.png"),
       p_residuals, width = 9, height = 7, dpi = 150)
ggsave(file.path(OUT_DIR, "standardized_residuals_cross_class.pdf"),
       p_std_resid, width = 8, height = 5)
ggsave(file.path(OUT_DIR, "standardized_residuals_cross_class.png"),
       p_std_resid, width = 8, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "slope_comparison_methods.pdf"),
       p_slopes, width = 8, height = 5)
ggsave(file.path(OUT_DIR, "slope_comparison_methods.png"),
       p_slopes, width = 8, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "reference_predictions_by_class.pdf"),
       p_refpred, width = 8, height = 5)
ggsave(file.path(OUT_DIR, "reference_predictions_by_class.png"),
       p_refpred, width = 8, height = 5, dpi = 150)

# ── FIGURES 5 & 6: BOXPLOTS BY ORDER AND FAMILY ───────────────────────────────
# Y-axis = std_residual (z-score) for direct cross-class comparison.

boxplot_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "none"
  )

order_counts <- yardstick_flagged %>%
  filter(!is.na(Order_final)) %>%
  count(Class_final, Order_final, name = "n_pairs") %>%
  filter(n_pairs >= 5)

order_medians <- yardstick_flagged %>%
  filter(!is.na(Order_final)) %>%
  semi_join(order_counts, by = c("Class_final", "Order_final")) %>%
  group_by(Class_final, Order_final) %>%
  summarize(med = median(std_residual, na.rm = TRUE), .groups = "drop") %>%
  arrange(Class_final, med) %>%
  mutate(order_rank = row_number())

plot_data_order <- yardstick_flagged %>%
  filter(!is.na(Order_final)) %>%
  semi_join(order_counts, by = c("Class_final", "Order_final")) %>%
  left_join(order_medians, by = c("Class_final", "Order_final")) %>%
  mutate(Order_label = sprintf(
    "%s (n=%d)", Order_final,
    order_counts$n_pairs[
      match(paste(Class_final, Order_final),
            paste(order_counts$Class_final, order_counts$Order_final))]))

p_order <- ggplot(plot_data_order,
                  aes(x = reorder(Order_label, order_rank),
                      y = std_residual,
                      fill  = Class_final,
                      color = Class_final)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted",
             color = "gray65", linewidth = 0.5) +
  geom_hline(yintercept = c(-2.576, 2.576), linetype = "dotted",
             color = "gray40", linewidth = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, linewidth = 0.5, width = 0.6) +
  geom_jitter(aes(color = Class_final),
              width = 0.2, height = 0, alpha = 0.45, size = 1.2, shape = 16) +
  geom_jitter(data = plot_data_order %>% filter(outlier_95),
              color = "black", fill = "yellow",
              width = 0.2, height = 0, alpha = 0.9, size = 2.0,
              shape = 21, stroke = 0.6) +
  geom_jitter(data = plot_data_order %>% filter(outlier_99),
              color = "black", fill = "red",
              width = 0.2, height = 0, alpha = 0.9, size = 2.4,
              shape = 21, stroke = 0.7) +
  scale_fill_manual(values  = class_colors) +
  scale_color_manual(values = class_colors) +
  facet_wrap(~ Class_final, scales = "free_x", ncol = 2) +
  labs(
    title    = "Time-corrected genomic divergence by Order",
    subtitle = "Std residual log(k2p) | dashed = expected | dotted = ±1.96 (95%) and ±2.576 (99%) | yellow = 95% outlier | red = 99% outlier",
    x = NULL, y = "Standardized residual (z-score)") +
  boxplot_theme

family_counts <- yardstick_flagged %>%
  filter(!is.na(Family_final)) %>%
  count(Class_final, Family_final, name = "n_pairs") %>%
  filter(n_pairs >= 5)

family_medians <- yardstick_flagged %>%
  filter(!is.na(Family_final)) %>%
  semi_join(family_counts, by = c("Class_final", "Family_final")) %>%
  group_by(Class_final, Family_final) %>%
  summarize(med = median(std_residual, na.rm = TRUE), .groups = "drop") %>%
  arrange(Class_final, med) %>%
  mutate(family_rank = row_number())

plot_data_family <- yardstick_flagged %>%
  filter(!is.na(Family_final)) %>%
  semi_join(family_counts, by = c("Class_final", "Family_final")) %>%
  left_join(family_medians, by = c("Class_final", "Family_final")) %>%
  mutate(Family_label = sprintf(
    "%s (n=%d)", Family_final,
    family_counts$n_pairs[
      match(paste(Class_final, Family_final),
            paste(family_counts$Class_final, family_counts$Family_final))]))

p_family <- ggplot(plot_data_family,
                   aes(x = reorder(Family_label, family_rank),
                       y = std_residual,
                       fill  = Class_final,
                       color = Class_final)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted",
             color = "gray65", linewidth = 0.5) +
  geom_hline(yintercept = c(-2.576, 2.576), linetype = "dotted",
             color = "gray40", linewidth = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, linewidth = 0.5, width = 0.6) +
  geom_jitter(aes(color = Class_final),
              width = 0.2, height = 0, alpha = 0.45, size = 1.0, shape = 16) +
  geom_jitter(data = plot_data_family %>% filter(outlier_95),
              color = "black", fill = "yellow",
              width = 0.2, height = 0, alpha = 0.9, size = 1.8,
              shape = 21, stroke = 0.6) +
  geom_jitter(data = plot_data_family %>% filter(outlier_99),
              color = "black", fill = "red",
              width = 0.2, height = 0, alpha = 0.9, size = 2.2,
              shape = 21, stroke = 0.7) +
  scale_fill_manual(values  = class_colors) +
  scale_color_manual(values = class_colors) +
  facet_wrap(~ Class_final, scales = "free_x", ncol = 2) +
  labs(
    title    = "Time-corrected genomic divergence by Family",
    subtitle = "Std residual log(k2p) | dashed = expected | dotted = ±1.96 (95%) and ±2.576 (99%) | yellow = 95% outlier | red = 99% outlier",
    x = NULL, y = "Standardized residual (z-score)") +
  boxplot_theme

ggsave(file.path(OUT_DIR, "divergence_by_order.pdf"),
       p_order, width = 14, height = 10)
ggsave(file.path(OUT_DIR, "divergence_by_order.png"),
       p_order, width = 14, height = 10, dpi = 200)

n_families <- n_distinct(plot_data_family$Family_final)
fig_width  <- max(14, ceiling(n_families * 0.35))
ggsave(file.path(OUT_DIR, "divergence_by_family.pdf"),
       p_family, width = fig_width, height = 11)
ggsave(file.path(OUT_DIR, "divergence_by_family.png"),
       p_family, width = fig_width, height = 11, dpi = 200)

cat(sprintf("\nOrder boxplot: %d orders across %d classes\n",
            n_distinct(plot_data_order$Order_final),
            n_distinct(plot_data_order$Class_final)))
cat(sprintf("Family boxplot: %d families across %d classes\n",
            n_distinct(plot_data_family$Family_final),
            n_distinct(plot_data_family$Class_final)))

# ── STEP 10: ANALYSIS DATA SUBSET (for PMM downstream) ────────────────────────
# Aves + Mammalia, 15 MYA cap — preserved from prior version so that
# downstream PMM code that loads this script's namespace still works.

analysis_data <- merged_data %>%
  filter(
    Class_final %in% c("Aves", "Mammalia"),
    !is.na(k2p),
    !is.na(timetree_div),
    k2p > 0,
    timetree_div > 0,
    timetree_div <= 15
  ) %>%
  mutate(
    Class  = factor(Class_final),
    Order  = factor(Order_final),
    Family = factor(Family_final),
    genus  = factor(genus_final)
  )

cat(sprintf("\nDownstream PMM analysis_data pairs: %d (Aves + Mammalia, ≤15 MYA)\n",
            nrow(analysis_data)))
print(table(analysis_data$Class_final))

# ── STEP 11: SHINY APP DATA EXPORT ────────────────────────────────────────────

SHINY_DIR <- "shiny_yardstick/data"
dir.create(SHINY_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Fitted model objects
saveRDS(class_fits,        file.path(SHINY_DIR, "class_fits.rds"))
saveRDS(class_fits_uncapped, file.path(SHINY_DIR, "class_fits_uncapped.rds"))
saveRDS(class_fits_common,  file.path(SHINY_DIR, "class_fits_common.rds"))
saveRDS(class_div_ranges,   file.path(SHINY_DIR, "class_div_ranges.rds"))
saveRDS(pooled_fit,         file.path(SHINY_DIR, "pooled_fit.rds"))
saveRDS(pooled_div_range,   file.path(SHINY_DIR, "pooled_div_range.rds"))
saveRDS(family_yardsticks,  file.path(SHINY_DIR, "family_yardsticks.rds"))

# 2. Pair-level public data (no genome accession IDs)
yardstick_public <- yardstick_flagged %>%
  dplyr::select(
    Class_final, Order_final, Family_final, genus_final,
    sp1, sp2, timetree_div, k2p, log_k2p, log_div_time,
    expected_log_k2p, pi_lower_95, pi_upper_95,
    pi_lower_99, pi_upper_99,
    residual, std_residual, p_value, q_fdr,
    outlier_sr_95, outlier_sr_95_low, outlier_sr_95_high,
    outlier_sr_99, outlier_sr_99_low, outlier_sr_99_high,
    outlier_95, outlier_95_low, outlier_95_high,
    outlier_99, outlier_99_low, outlier_99_high,
    outlier_fdr,
    pooled_expected_log, pooled_pi_lower, pooled_pi_upper,
    pooled_residual, pooled_std_residual,
    pooled_outlier, pooled_outlier_low, pooled_outlier_high
  )
saveRDS(yardstick_public, file.path(SHINY_DIR, "yardstick_data.rds"))

cat(sprintf("\nShiny app data saved to: %s/\n", SHINY_DIR))
cat(sprintf("\nYardstick analysis complete. All outputs in: %s/\n", OUT_DIR))

# ── SUMMARY ───────────────────────────────────────────────────────────────────

cat("\n=========================================================================\n")
cat("YARDSTICK ANALYSIS SUMMARY\n")
cat("=========================================================================\n")
cat(sprintf("Primary fits:      genus-capped (≤ %d pairs/genus)\n", GENUS_CAP))
cat(sprintf("Sensitivities:     uncapped OLS, common-window (%g–%g MYA), LMM\n",
            COMMON_WINDOW[1], COMMON_WINDOW[2]))
cat(sprintf("Diagnostics:       Shapiro–Wilk, Breusch–Pagan, segmented regression\n"))
cat(sprintf("Outlier tiers:     std resid |z|>1.96/2.576 (primary), 95%%/99%% PI, BH-FDR (q < %.2f)\n", FDR_ALPHA))
cat(sprintf("Outlier cross-chk: DHARMa testOutliers() on class OLS + genus LMM\n"))
cat(sprintf("Family yardstick:  two-tier (r² ≥ %.3f [mode=%s] OLS vs. LMM intercept)\n",
            r2_cutoff_used, R2_CUTOFF_MODE))
cat(sprintf("Cross-class:       reference predictions at %s MYA\n",
            paste(REF_TIMES_MYA, collapse = ", ")))
cat("=========================================================================\n")

# =============================================================================
# ADD-ON: REPORT THE ALL-CLASS (POOLED) YARDSTICK + SURFACE AMPHIBIA
# =============================================================================
# Why amphibians "disappear": fit_class_regression() returns NULL when a class
# has < MIN_CLASS_PAIRS (20) pairs after the >=80% alignment filter. Amphibian
# genomes align poorly, so Amphibia usually falls below that threshold and gets
# no class-specific fit. It therefore drops out of:
#     class_fits, class_yardstick_table, diagnostics, reference_predictions,
#     and the class figures.
# The pooled (all-class) fit IS computed (pooled_fit, ~line 457) but is only
# printed with cat() -- it never enters the comparison table, the prediction
# table, or any figure. This add-on promotes the pooled model to a first-class
# citizen and reports Amphibia through it.
#
# Run AFTER the existing yardstick_analysis.R has built: pooled_fit,
# pooled_sigma, slope_comparison, reference_predictions, yardstick_flagged,
# class_colors, REF_TIMES_MYA, OUT_DIR, yardstick_data(_capped).
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1) Add the pooled model to the slope/intercept comparison table
# -----------------------------------------------------------------------------
pooled_row <- {
  s <- summary(pooled_fit)
  tibble(
    class       = "All classes (pooled)",
    method      = "pooled_capped",
    n           = length(residuals(pooled_fit)),
    intercept   = unname(coef(pooled_fit)[1]),
    slope       = unname(coef(pooled_fit)[2]),
    slope_se    = unname(coef(s)[2, "Std. Error"]),
    r_squared   = s$r.squared,
    residual_sd = s$sigma
  )
}

slope_comparison_with_pooled <- bind_rows(slope_comparison, pooled_row) %>%
  arrange(class, method)

cat("\n=== SLOPE / INTERCEPT COMPARISON (incl. pooled all-class) ===\n")
print(slope_comparison_with_pooled, n = 40)
write_csv(slope_comparison_with_pooled,
          file.path(OUT_DIR, "slope_comparison_methods_with_pooled.csv"))

# -----------------------------------------------------------------------------
# 2) Pooled reference predictions at the standard divergence times
#    (so the all-class model appears alongside the per-class columns)
# -----------------------------------------------------------------------------
pooled_ref <- tibble(div_time_mya = REF_TIMES_MYA) %>%
  mutate(
    .pred = map(div_time_mya, ~ {
      nd <- tibble(log_div_time = log(.x))
      pr <- predict(pooled_fit, newdata = nd, interval = "prediction", level = 0.95)
      tibble(expected_k2p = exp(pr[, "fit"]),
             lower = exp(pr[, "lwr"]), upper = exp(pr[, "upr"]))
    })
  ) %>%
  unnest(.pred) %>%
  mutate(class = "All classes (pooled)",
         k2p_pretty = sprintf("%.4f [%.4f - %.4f]", expected_k2p, lower, upper))

cat("\n=== POOLED ALL-CLASS REFERENCE PREDICTIONS ===\n")
print(pooled_ref %>% dplyr::select(div_time_mya, k2p_pretty))
write_csv(pooled_ref, file.path(OUT_DIR, "reference_predictions_pooled.csv"))

# -----------------------------------------------------------------------------
# 3) Amphibian results, reported THROUGH the pooled model
#    yardstick_flagged already carries pooled_* columns for every pair; we just
#    subset to Amphibia and summarize, since they have no class-specific fit.
# -----------------------------------------------------------------------------
amphibia_pooled <- yardstick_flagged %>%
  filter(Class_final == "Amphibia") %>%
  mutate(
    residual_pooled = pooled_residual,
    std_resid_pooled = pooled_std_residual
  )

amphibia_summary <- amphibia_pooled %>%
  summarize(
    class            = "Amphibia",
    model            = "pooled all-class",
    n_pairs          = n(),
    n_genera         = n_distinct(genus_final),
    n_families       = n_distinct(Family_final),
    median_divtime   = median(timetree_div),
    range_divtime    = sprintf("%.2f - %.2f My", min(timetree_div), max(timetree_div)),
    resid_mean       = mean(pooled_residual,    na.rm = TRUE),
    resid_sd         = sd(pooled_residual,      na.rm = TRUE),
    n_outlier_pooled = sum(pooled_outlier,      na.rm = TRUE)
  )

cat("\n=== AMPHIBIA (via pooled all-class model) ===\n")
print(amphibia_summary)
write_csv(amphibia_summary, file.path(OUT_DIR, "amphibia_pooled_summary.csv"))

# Amphibian outlier pairs flagged against the pooled PI
write_csv(
  amphibia_pooled %>%
    filter(pooled_outlier) %>%
    dplyr::select(Class_final, Family_final, genus_final, sp1, sp2,
                  timetree_div, k2p,
                  pooled_expected_log, pooled_pi_lower, pooled_pi_upper,
                  pooled_residual, pooled_std_residual,
                  pooled_outlier_low, pooled_outlier_high),
  file.path(OUT_DIR, "amphibia_outlier_pairs_pooled.csv")
)

# -----------------------------------------------------------------------------
# 4) Figure: divergence curve for the POOLED model with all four classes' raw
#    points overlaid (amphibians included). Complements the per-class figure.
# -----------------------------------------------------------------------------
pooled_range <- range(yardstick_data_capped$timetree_div)
time_seq <- exp(seq(log(pooled_range[1]), log(pooled_range[2]), length.out = 200))
pooled_ribbon <- {
  pr <- predict(pooled_fit, newdata = tibble(log_div_time = log(time_seq)),
                interval = "prediction")
  tibble(div_time = time_seq,
         fit = exp(pr[, "fit"]), lower = exp(pr[, "lwr"]), upper = exp(pr[, "upr"]))
}

p_pooled <- ggplot() +
  geom_ribbon(data = pooled_ribbon,
              aes(x = div_time, ymin = lower, ymax = upper),
              fill = "gray70", alpha = 0.3) +
  geom_line(data = pooled_ribbon, aes(x = div_time, y = fit),
            color = "black", linewidth = 1.1) +
  geom_point(data = yardstick_data,
             aes(x = timetree_div, y = k2p, color = Class_final),
             alpha = 0.35, size = 0.9) +
  scale_color_manual(values = class_colors) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title    = "Pooled all-class genomic yardstick",
    subtitle = "Black line/grey ribbon = pooled fit + 95% PI; points colored by Class (incl. Amphibia)",
    x = "Divergence time (My, log scale)", y = "k2p divergence (log scale)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "pooled_allclass_yardstick.png"),
       p_pooled, width = 9, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "pooled_allclass_yardstick.pdf"),
       p_pooled, width = 9, height = 6)
cat("\nSaved pooled all-class yardstick figure and tables.\n")

#### plot a better boxplot figure
# ============================================================================
# Publication figure: time-corrected genomic divergence by Order
#   * horizontal boxplots, ONE CLASS PER ROW (panel height ∝ #orders, so every
#     box is the same size — fixes the unequal-width problem)
#   * Birds = blue, Mammals = red, Reptiles = green
#   * flagged 95% / 99% prediction-interval outlier pairs highlighted
#
# INPUT: the `yardstick_flagged` data frame from the yardstick pipeline. It must
# contain the columns: Class_final, Order_final, std_residual, outlier_95,
# outlier_99. If you saved it (saveRDS(yardstick_flagged, "yardstick_flagged.rds"))
# set INPUT_RDS below; otherwise the script uses the in-memory object.
# ============================================================================

library(tidyverse)

# ---- settings ---------------------------------------------------------------
INPUT_RDS <- "yardstick_flagged.rds"   # used only if `yardstick_flagged` isn't in memory
MIN_N     <- 5                          # minimum pairs per Order to include
OUT_STEM  <- "divergence_by_order_pub"  # output filename stem

dat <- if (exists("yardstick_flagged")) yardstick_flagged else readRDS(INPUT_RDS)

# Aves→Birds etc.; blue / red / green
class_map  <- c(Aves = "Birds", Mammalia = "Mammals", Reptilia = "Reptiles")
class_cols <- c(Birds = "#2C6FB3", Mammals = "#C0392B", Reptiles = "#2E8B57")

# NOTE: dplyr verbs are namespaced (dplyr::) so this works even inside the
# yardstick script, where other attached packages mask select()/filter().

# ---- prepare data -----------------------------------------------------------
d <- dat %>%
  dplyr::filter(Class_final %in% names(class_map),
                !is.na(Order_final), !is.na(std_residual)) %>%
  dplyr::mutate(Class = factor(class_map[as.character(Class_final)],
                               levels = c("Birds", "Mammals", "Reptiles")))

# keep Orders with enough pairs, attach n
counts <- d %>% dplyr::count(Class, Order_final, name = "n") %>%
  dplyr::filter(n >= MIN_N)
d <- d %>% dplyr::inner_join(counts, by = c("Class", "Order_final"))

# order Orders within each class by median residual; label with sample size
ord <- d %>%
  dplyr::group_by(Class, Order_final, n) %>%
  dplyr::summarise(med = median(std_residual), .groups = "drop") %>%
  dplyr::arrange(Class, med) %>%
  dplyr::mutate(label = sprintf("%s (n = %d)", Order_final, n))
d <- d %>%
  dplyr::left_join(dplyr::select(ord, Class, Order_final, label),
                   by = c("Class", "Order_final")) %>%
  dplyr::mutate(label = factor(label, levels = ord$label))  # class-blocked, sorted by median

# outlier points, tiered (99% takes precedence over 95%)
outliers <- d %>%
  dplyr::filter(outlier_95 | outlier_99) %>%
  dplyr::mutate(tier = factor(ifelse(outlier_99, "Beyond 99% PI", "Beyond 95% PI"),
                              levels = c("Beyond 95% PI", "Beyond 99% PI")))

# ---- plot -------------------------------------------------------------------
p <- ggplot(d, aes(x = std_residual, y = label)) +
  # reference lines: expected (0) and the ±1.96 / ±2.576 standardized bands
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
  geom_vline(xintercept = c(-1.96, 1.96), colour = "grey70",
             linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = c(-2.576, 2.576), colour = "grey80",
             linetype = "dotted", linewidth = 0.4) +
  # boxplots, filled by class (outline fixed so the colour scale is free for points)
  geom_boxplot(aes(fill = Class), colour = "grey25",
               width = 0.62, alpha = 0.45, linewidth = 0.45,
               outlier.shape = NA) +
  # highlighted outlier pairs
  geom_point(data = outliers,
             aes(colour = tier, shape = tier, size = tier),
             position = position_jitter(height = 0.18, width = 0),
             stroke = 0.6, alpha = 0.95) +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_colour_manual(values = c("Beyond 95% PI" = "grey25",
                                 "Beyond 99% PI" = "#D35400"), name = NULL) +
  scale_shape_manual(values = c("Beyond 95% PI" = 1, "Beyond 99% PI" = 18), name = NULL) +
  scale_size_manual(values  = c("Beyond 95% PI" = 1.7, "Beyond 99% PI" = 2.9), name = NULL) +
  facet_grid(Class ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Time-corrected genomic divergence by Order",
    subtitle = "Standardized residual from each Class's divergence–time regression; boxes ordered by median",
    x        = "Time-corrected divergence (standardized residual of log k2p)",
    y        = NULL,
    caption  = "Dashed = ±1.96, dotted = ±2.576. Points: pairs outside the 95% / 99% prediction interval."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 15),
    plot.subtitle   = element_text(size = 10.5, colour = "grey30", margin = margin(b = 6)),
    plot.caption    = element_text(size = 9, colour = "grey40", hjust = 0),
    axis.text.y     = element_text(size = 10, colour = "grey10"),
    axis.text.x     = element_text(size = 11),
    axis.title.x    = element_text(size = 12, margin = margin(t = 8)),
    strip.text.y    = element_text(face = "bold", size = 12, angle = 0),
    strip.background = element_rect(fill = "grey93", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.spacing   = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.margin   = margin(t = -2),
    plot.margin     = margin(10, 16, 8, 8)
  )

# ---- save (height scales with the number of Orders) -------------------------
n_orders <- nlevels(droplevels(d$label))
fig_h    <- max(5, 0.42 * n_orders + 1.6)

# PDF: the base `pdf` device is always available (no Cairo/X11 needed).
pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(paste0(OUT_STEM, ".pdf"), p, width = 8.5, height = fig_h, device = pdf_dev)
cat(sprintf("Saved %s.pdf (%d orders, %.1f in tall)\n", OUT_STEM, n_orders, fig_h))

# PNG: needs a raster device. On headless HPC that means ragg or Cairo; if
# neither exists, don't halt — the PDF is the publication file, and you can
# rasterize it later (e.g. `pdftoppm -r 400 -png <pdf>`).
png_ok <- tryCatch({
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(paste0(OUT_STEM, ".png"), p, width = 8.5, height = fig_h, dpi = 400,
           device = ragg::agg_png)
  } else if (capabilities("cairo")) {
    ggsave(paste0(OUT_STEM, ".png"), p, width = 8.5, height = fig_h, dpi = 400,
           type = "cairo")
  } else {
    ggsave(paste0(OUT_STEM, ".png"), p, width = 8.5, height = fig_h, dpi = 400)
  }
  TRUE
}, error = function(e) { cat(sprintf("  PNG skipped (%s). PDF is fine.\n", conditionMessage(e))); FALSE })
if (png_ok) cat(sprintf("Saved %s.png\n", OUT_STEM))

saveRDS(yardstick_flagged, "yardstick_flagged.rds")

