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
library(lme4)
library(lmtest)
library(segmented)
library(DHARMa)

# mito yardstick


MIN_QUERY_ALIGN     <- 80
MIN_REF_ALIGN       <- 80
MIN_CLASS_PAIRS     <- 20
GENUS_CAP           <- 25
COMMON_WINDOW       <- c(5, 15)
PI_LEVELS           <- c(0.95, 0.99)
FDR_ALPHA           <- 0.05
REF_TIMES_MYA       <- c(1, 2, 5, 10, 15)
SEED                <- 12345

OUT_DIR <- "mito_yardstick"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "diagnostics"), showWarnings = FALSE)

# load data

setwd("/scratch/gautschi/allen715/GD_mito")
traits     <- read_excel("/scratch/gautschi/allen715/GD_models/final_dataset/Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
divergence <- read_csv("/scratch/gautschi/allen715/GD_models/final_dataset/Complete_divergence_with_mito_05-14-26.csv",
                       guess_max = Inf, show_col_types = FALSE) %>%
  # force genome-accession keys to character (avoids trait-join type mismatch)
  dplyr::mutate(dplyr::across(dplyr::any_of(c("genome1", "genome2")), as.character))

colnames(traits)     <- make.names(colnames(traits))
colnames(divergence) <- make.names(colnames(divergence))

# mito: response = coding 13-pcg k2p (aliased to k2p); QC-flagged pairs dropped; 80% coverage kept for parity
divergence <- divergence %>%
  filter(is.na(mito_k2p_flag) | mito_k2p_flag == "") %>%
  mutate(
    k2p                  = suppressWarnings(as.numeric(mito_k2p)),
    mito_query_align_pct = suppressWarnings(as.numeric(mito_query_align_pct)),
    mito_ref_align_pct   = suppressWarnings(as.numeric(mito_ref_align_pct))
  )

cat(sprintf("Loaded %d species with traits\n",        nrow(traits)))
cat(sprintf("Loaded %d species pairs with divergence data\n", nrow(divergence)))

### parse genome metrics

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

### exclude flagged species

excluded_species <- traits %>%
  filter(Exclude. %in% c("yes", "Yes")) %>%
  pull(accession)

flag_domestic_pair <- function(div) {
  hits <- div %>%
    filter(
      (str_detect(species1, "frontalis") & str_detect(species2, "gaurus")) |
      (str_detect(species1, "gaurus") & str_detect(species2, "frontalis"))
    )
  if (nrow(hits) > 0) {

  }
}
flag_domestic_pair(divergence)


divergence_clean <- divergence %>%
  filter(!genome1 %in% excluded_species,
         !genome2 %in% excluded_species)


### alignment quality filter

divergence_clean <- divergence_clean %>%
  filter(
    !is.na(mito_query_align_pct),
    !is.na(mito_ref_align_pct),
    mito_query_align_pct >= MIN_QUERY_ALIGN,
    mito_ref_align_pct   >= MIN_REF_ALIGN
  )

cat(sprintf("Pairs after alignment filters (>=%d%% both genomes): %d\n",
            MIN_QUERY_ALIGN, nrow(divergence_clean)))

### name normalisation

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

### randomize sp assignment

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

### merge traits to pairs

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

# pair-level trait averages and differences
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

# yardstick dataset

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

# genus cap

set.seed(SEED)
yardstick_data_capped <- yardstick_data %>%
  group_by(Class_final, genus_final) %>%
  slice_sample(n = GENUS_CAP, replace = FALSE) %>%
  ungroup()

write_csv(yardstick_data, file.path(OUT_DIR, "yardstick_pairs_all.csv"))

genus_cap_summary <- yardstick_data %>%
  group_by(Class_final, genus_final) %>%
  summarize(n_before = n(), .groups = "drop") %>%
  mutate(n_after = pmin(n_before, GENUS_CAP),
         n_dropped = n_before - n_after) %>%
  filter(n_dropped > 0) %>%
  arrange(desc(n_dropped))

if (nrow(genus_cap_summary) > 0) {
  cat("Top genera capped (n_before -> n_after):\n")
  print(genus_cap_summary %>% head(10), n = 10)
}

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

# class level regressions

fit_class_regression <- function(df, class_name, label = "") {
  df_class <- df %>% filter(Class_final == class_name)

  if (nrow(df_class) < MIN_CLASS_PAIRS) {
    cat(sprintf("  %s [%s]: insufficient pairs (%d) — skipping\n",
                class_name, label, nrow(df_class)))
    return(NULL)
  }

  fit <- lm(log_k2p ~ log_div_time, data = df_class)

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


class_fits <- map(
  levels(yardstick_data_capped$Class_final),
  ~ fit_class_regression(yardstick_data_capped, .x, label = "primary, capped")
) %>%
  setNames(levels(yardstick_data_capped$Class_final))

# per-class divergence-time ranges
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

# sensitivity fits

# uncapped ols
class_fits_uncapped <- map(
  levels(yardstick_data$Class_final),
  ~ fit_class_regression(yardstick_data, .x, label = "uncapped")
) %>%
  setNames(levels(yardstick_data$Class_final))

# common divergence-time window (5–15 MYA)
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

# LMM with genus random intercept
cat("\n=== SENSITIVITY 3: LMM (RANDOM INTERCEPT BY GENUS) ===")
fit_class_lmm <- function(df, class_name) {
  df_class <- df %>% filter(Class_final == class_name)
  if (nrow(df_class) < MIN_CLASS_PAIRS) return(NULL)
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

# pooled all-class fit
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

# comparison table

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
    r_squared   = NA_real_,
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

print(slope_comparison, n = 30)
write_csv(slope_comparison, file.path(OUT_DIR, "slope_comparison_methods.csv"))

# residual diagnostics


diagnose_class <- function(fit, df_class, class_name) {
  if (is.null(fit)) return(NULL)
  resid_vec <- residuals(fit)
  fitted_vec <- fitted(fit)

  # normality
  sw <- if (length(resid_vec) >= 3 && length(resid_vec) <= 5000) {
    shapiro.test(resid_vec)
  } else NULL

  # heteroscedasticity (Breusch-Pagan)
  bp <- tryCatch(bptest(fit), error = function(e) NULL)

  # time-dependence: segmented regression
  seg <- tryCatch(
    segmented(fit, seg.Z = ~ log_div_time),
    error = function(e) NULL
  )
  seg_breakpoint_mya <- if (!is.null(seg) && !is.null(seg$psi)) {
    exp(seg$psi[1, "Est."])
  } else NA_real_
  # davies test for non-constant slope
  seg_pvalue <- tryCatch(
    davies.test(fit, seg.Z = ~ log_div_time)$p.value,
    error = function(e) NA_real_
  )

  # binned residual sd
  bin_summary <- tibble(fitted = fitted_vec, residual = resid_vec) %>%
    mutate(bin = cut(fitted, breaks = quantile(fitted, probs = seq(0, 1, 0.2),
                                                na.rm = TRUE),
                    include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarize(n = n(), bin_sd = sd(residual), .groups = "drop")

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

# qq plots
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

# compute yardsticks

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
print(yardstick_table)
write_csv(yardstick_table, file.path(OUT_DIR, "class_yardstick_table.csv"))

# predict yardstick for new pair

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

# reference predictions


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

ref_wide <- reference_predictions %>%
  dplyr::select(div_time_mya, class, k2p_pretty) %>%
  pivot_wider(names_from = class, values_from = k2p_pretty)

cat("Predicted k2p [95% PI] by class at standard divergence times:\n")
print(ref_wide)

write_csv(reference_predictions, file.path(OUT_DIR, "reference_predictions_by_class.csv"))
write_csv(ref_wide, file.path(OUT_DIR, "reference_predictions_by_class_wide.csv"))

# outlier flagging

flag_outliers_extended <- function(df, class_fits) {
  df %>%
    group_by(Class_final) %>%
    group_modify(~ {
      cn  <- as.character(unique(.y$Class_final))
      fit <- class_fits[[cn]]
      if (is.null(fit)) return(.x)

      sigma  <- summary(fit)$sigma
      df_fit <- length(residuals(fit)) - 2

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

# pooled-model flags for Amphibia
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

print(outlier_summary)
write_csv(outlier_summary, file.path(OUT_DIR, "outlier_summary_by_class.csv"))

# flagged pairs csv (95%)
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

# 99% pairs
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

# standardized-residual pairs (|z| > 1.96)
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

# dharma cross-check of outliers


set.seed(SEED)

run_dharma_outlier_test <- function(fit, class_name, model_label) {
  if (is.null(fit)) return(NULL)
  res <- tryCatch({
    sim <- simulateResiduals(fittedModel = fit, n = 250, plot = FALSE)
    ot  <- testOutliers(sim, type = "binomial", plot = FALSE)
    du  <- testUniformity(sim, plot = FALSE)  # KS test on residual uniformity
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
  print(dharma_outlier_results)
  write_csv(dharma_outlier_results,
            file.path(OUT_DIR, "diagnostics", "dharma_outlier_tests.csv"))
}

# clade-level outlier summary

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

print(genus_outlier_summary, n = 50)
write_csv(genus_outlier_summary,
          file.path(OUT_DIR, "genus_outlier_summary.csv"))

# figures

class_colors <- c(
  "Aves"     = "#2980B9",
  "Mammalia" = "#C0392B",
  "Reptilia" = "#27AE60",
  "Amphibia" = "#9e22e6"
)

# class-level yardstick curves
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

# raw residual densities per class
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

# standardized residuals
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

# slope comparison across methods
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

# reference-time predictions (johns & avise comparison)
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

### pooled all class yardstick

library(tidyverse)

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

print(slope_comparison_with_pooled, n = 40)
write_csv(slope_comparison_with_pooled,
          file.path(OUT_DIR, "slope_comparison_methods_with_pooled.csv"))

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

print(pooled_ref %>% dplyr::select(div_time_mya, k2p_pretty))
write_csv(pooled_ref, file.path(OUT_DIR, "reference_predictions_pooled.csv"))

### amphibian results
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
    median_divtime   = median(timetree_div),
    range_divtime    = sprintf("%.2f - %.2f My", min(timetree_div), max(timetree_div)),
    resid_mean       = mean(pooled_residual,    na.rm = TRUE),
    resid_sd         = sd(pooled_residual,      na.rm = TRUE),
    n_outlier_pooled = sum(pooled_outlier,      na.rm = TRUE)
  )

print(amphibia_summary)
write_csv(amphibia_summary, file.path(OUT_DIR, "amphibia_pooled_summary.csv"))

write_csv(
  amphibia_pooled %>%
    filter(pooled_outlier) %>%
    dplyr::select(Class_final, genus_final, sp1, sp2,
                  timetree_div, k2p,
                  pooled_expected_log, pooled_pi_lower, pooled_pi_upper,
                  pooled_residual, pooled_std_residual,
                  pooled_outlier_low, pooled_outlier_high),
  file.path(OUT_DIR, "amphibia_outlier_pairs_pooled.csv")
)

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
