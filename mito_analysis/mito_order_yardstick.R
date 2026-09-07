### mitochondrial order-level yardstick

library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(readxl)
library(readr)
library(tidyverse)
library(lme4) 

MIN_QUERY_ALIGN  <- 80
MIN_REF_ALIGN    <- 80
GENUS_CAP        <- 25  
MIN_ORDER_PAIRS  <- 15  
MIN_ORDER_TIMES  <- 3  
MIN_CLASS_PAIRS  <- 20  
USE_GENUS_CAP    <- TRUE  
R2_THRESHOLD     <- 0.30  
R2_CUTOFF_MODE   <- "quantile"  # "fixed" will use threshold
R2_QUANTILE_TARGET <- 0.50
FDR_ALPHA        <- 0.05
REF_TIMES_MYA    <- c(1, 2, 5, 10, 15)
SEED             <- 12345

OUT_DIR <- "mito_order_yardstick"
dir.create(file.path(OUT_DIR, "diagnostics"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"),     showWarnings = FALSE, recursive = TRUE)

setwd("/scratch/gautschi/allen715/GD_mito")

# load data
traits     <- read_excel("/scratch/gautschi/allen715/GD_models/final_dataset/Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
divergence <- read_csv("/scratch/gautschi/allen715/GD_models/final_dataset/Complete_divergence_with_mito_05-14-26.csv",
                       guess_max = Inf, show_col_types = FALSE) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(c("genome1", "genome2")), as.character))

colnames(traits)     <- make.names(colnames(traits))
colnames(divergence) <- make.names(colnames(divergence))

divergence <- divergence %>%
  filter(is.na(mito_k2p_flag) | mito_k2p_flag == "") %>%  
  mutate(
    k2p                  = suppressWarnings(as.numeric(mito_k2p)),
    mito_query_align_pct = suppressWarnings(as.numeric(mito_query_align_pct)),
    mito_ref_align_pct   = suppressWarnings(as.numeric(mito_ref_align_pct))
  )

# exclude flagged species
excluded_species <- traits %>%
  filter(Exclude. %in% c("yes", "Yes")) %>% pull(accession)

divergence_clean <- divergence %>%
  filter(!genome1 %in% excluded_species, !genome2 %in% excluded_species)

# alignment quality filter
divergence_clean <- divergence_clean %>%
  filter(!is.na(mito_query_align_pct), !is.na(mito_ref_align_pct),
         mito_query_align_pct >= MIN_QUERY_ALIGN,
         mito_ref_align_pct   >= MIN_REF_ALIGN)

# name normalisation
normalize_species <- function(x) {
  x %>% str_trim() %>% str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>% str_replace_all("[^A-Za-z ]", "") %>%
    str_to_lower() %>% str_replace_all(" ", "_")
}
divergence_clean <- divergence_clean %>%
  mutate(species1_norm = normalize_species(species1),
         species2_norm = normalize_species(species2))

# randomize sp assignment
set.seed(SEED)
divergence_clean <- divergence_clean %>%
  rowwise() %>%
  mutate(swap       = sample(c(TRUE, FALSE), 1),
         sp1        = ifelse(swap, species2_norm, species1_norm),
         sp2        = ifelse(swap, species1_norm, species2_norm),
         genome_sp1 = ifelse(swap, genome2, genome1),
         genome_sp2 = ifelse(swap, genome1, genome2)) %>%
  ungroup() %>% dplyr::select(-swap)

# merge taxonomy to pairs 
taxo <- traits %>%
  dplyr::select(accession, Class, Order, Family, genus) %>%
  distinct(accession, .keep_all = TRUE)

merged_data <- divergence_clean %>%
  left_join(taxo %>% rename_with(~ paste0(.x, "_sp1")),
            by = c("genome_sp1" = "accession_sp1")) %>%
  left_join(taxo %>% rename_with(~ paste0(.x, "_sp2")),
            by = c("genome_sp2" = "accession_sp2")) %>%
  mutate(Class_final  = coalesce(Class_sp1,  Class_sp2),
         Order_final  = coalesce(Order_sp1,  Order_sp2),
         Family_final = coalesce(Family_sp1, Family_sp2),
         genus_final  = str_to_lower(coalesce(genus_sp1, genus_sp2)))

# build yardstick dataset 
yardstick_data <- merged_data %>%
  filter(Class_final %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"),
         !is.na(k2p), !is.na(timetree_div), k2p > 0, timetree_div > 0) %>%
  mutate(log_k2p      = log(k2p),
         log_div_time = log(timetree_div),
         Class_final  = factor(Class_final,
                               levels = c("Aves", "Mammalia", "Reptilia", "Amphibia")))

# genus cap 
set.seed(SEED)
yardstick_data_capped <- yardstick_data %>%
  group_by(Class_final, genus_final) %>%
  slice_sample(n = GENUS_CAP, replace = FALSE) %>% ungroup()

FIT_DATA <- if (USE_GENUS_CAP) yardstick_data_capped else yardstick_data

cat(sprintf("Pairs: %d total | fitting on %s data (%d pairs)\n",
            nrow(yardstick_data), if (USE_GENUS_CAP) "genus-capped" else "full",
            nrow(FIT_DATA)))
print(yardstick_data %>% count(Class_final, name = "n_pairs"))

# order-level ols fits 
fit_order_regression <- function(df, class_name, order_name) {
  d <- df %>% filter(Class_final == class_name, Order_final == order_name)
  if (nrow(d) < MIN_ORDER_PAIRS)               return(NULL)
  if (n_distinct(d$timetree_div) < MIN_ORDER_TIMES) return(NULL)
  fit <- tryCatch(lm(log_k2p ~ log_div_time, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  tibble(
    class = class_name, order = order_name,
    n_pairs = nrow(d), n_families = n_distinct(d$Family_final),
    n_genera = n_distinct(d$genus_final),
    intercept = coef(fit)[1], slope = coef(fit)[2],
    slope_se = s$coefficients[2, "Std. Error"],
    r_squared = s$r.squared, residual_sd = s$sigma,
    div_time_min = min(d$timetree_div), div_time_max = max(d$timetree_div)
  )
}

order_raw <- FIT_DATA %>%
  distinct(Class_final, Order_final) %>% filter(!is.na(Order_final)) %>%
  pmap_dfr(function(Class_final, Order_final)
    fit_order_regression(FIT_DATA, Class_final, Order_final))

cat(sprintf("\nOrders with a fit (n >= %d, >= %d distinct times): %d\n",
            MIN_ORDER_PAIRS, MIN_ORDER_TIMES, nrow(order_raw)))

# r2 cutooff
r2_vec <- order_raw$r_squared[is.finite(order_raw$r_squared)]
r2_quantile_table <- tibble(
  quantile  = c(0.10, 0.25, 0.50, 0.70, 0.75, 0.90),
  r_squared = quantile(r2_vec, c(0.10, 0.25, 0.50, 0.70, 0.75, 0.90), na.rm = TRUE))
r2_threshold_pctile <- 100 * ecdf(r2_vec)(R2_THRESHOLD)
r2_quantile_cutoff  <- as.numeric(quantile(r2_vec, R2_QUANTILE_TARGET, na.rm = TRUE))
r2_cutoff_used <- if (identical(R2_CUTOFF_MODE, "quantile")) r2_quantile_cutoff else R2_THRESHOLD

print(r2_quantile_table)

write_csv(r2_quantile_table, file.path(OUT_DIR, "diagnostics", "order_r2_quantiles.csv"))

order_ols <- order_raw %>%
  mutate(tier = if_else(r_squared >= r2_cutoff_used, "time_informative", "intercept_only"))
cat(sprintf("  time-informative: %d | intercept-only: %d\n",
            sum(order_ols$tier == "time_informative"),
            sum(order_ols$tier == "intercept_only")))

# class lmm with order random intercept 
fit_class_order_lmm <- function(df, class_name) {
  d <- df %>% filter(Class_final == class_name)
  if (nrow(d) < MIN_CLASS_PAIRS)          return(NULL)
  if (n_distinct(d$Order_final) < 3)      return(NULL)
  tryCatch(lmer(log_k2p ~ log_div_time + (1 | Order_final), data = d, REML = TRUE),
           error = function(e) NULL)
}
class_order_lmms <- map(levels(yardstick_data$Class_final),
                        ~ fit_class_order_lmm(yardstick_data, .x)) %>%
  setNames(levels(yardstick_data$Class_final))

extract_order_intercepts <- function(lmm, class_name) {
  if (is.null(lmm)) return(NULL)
  fe <- fixef(lmm); ints <- ranef(lmm)$Order_final
  tibble(class = class_name, order = rownames(ints),
         class_intercept = unname(fe[1]), class_slope = unname(fe[2]),
         order_intercept = unname(fe[1]) + ints[, "(Intercept)"],
         lmm_residual_sd = sigma(lmm))
}
order_lmm_intercepts <- map2_dfr(class_order_lmms, names(class_order_lmms),
                                 extract_order_intercepts)

# combine into final order yardstick 
order_yardsticks <- order_ols %>%
  left_join(order_lmm_intercepts, by = c("class", "order")) %>%
  mutate(
    final_intercept = if_else(tier == "time_informative", intercept, order_intercept),
    final_slope     = if_else(tier == "time_informative", slope,     class_slope),
    final_resid_sd  = if_else(tier == "time_informative", residual_sd, lmm_residual_sd)
  )

cat("\nOrder yardstick table:\n")
print(order_yardsticks %>% arrange(class, desc(n_pairs)) %>%
        dplyr::select(class, order, n_pairs, n_families, tier,
                      slope, r_squared, final_intercept, final_slope, final_resid_sd),
      n = 60)
write_csv(order_yardsticks, file.path(OUT_DIR, "order_yardstick_table.csv"))

# reference predictions per order
order_reference_predictions <- order_yardsticks %>%
  dplyr::select(class, order, final_intercept, final_slope, final_resid_sd,
                div_time_min, div_time_max) %>%
  tidyr::crossing(div_time_mya = REF_TIMES_MYA) %>%
  mutate(
    expected_log = final_intercept + final_slope * log(div_time_mya),
    lower_log    = expected_log - 1.96 * final_resid_sd,
    upper_log    = expected_log + 1.96 * final_resid_sd,
    expected_k2p = exp(expected_log), lower = exp(lower_log), upper = exp(upper_log),
    extrapolation = div_time_mya < div_time_min | div_time_mya > div_time_max
  ) %>%
  dplyr::select(class, order, div_time_mya, expected_k2p, lower, upper,
                expected_log, lower_log, upper_log, extrapolation)
write_csv(order_reference_predictions,
          file.path(OUT_DIR, "order_reference_predictions.csv"))

# outlier flagging
# class yardsticks fit on the genus-capped data (if you want to see the class outliers in the non-time-informative orders)
CLASS_YARDSTICK_CSV <- ""  
if (nzchar(CLASS_YARDSTICK_CSV) && file.exists(CLASS_YARDSTICK_CSV)) {
  class_yardsticks <- read_csv(CLASS_YARDSTICK_CSV, show_col_types = FALSE) %>%
    transmute(Class_final = class, class_intercept = intercept,
              class_slope = slope, class_resid_sd = residual_sd)
} else {
  fit_class_yardstick <- function(df, cn) {
    d <- df %>% filter(Class_final == cn)
    if (nrow(d) < MIN_CLASS_PAIRS) return(NULL)
    f <- lm(log_k2p ~ log_div_time, data = d); s <- summary(f)
    tibble(Class_final = cn, class_intercept = coef(f)[1],
           class_slope = coef(f)[2], class_resid_sd = s$sigma)
  }
  class_yardsticks <- map_dfr(levels(yardstick_data$Class_final),
                              ~ fit_class_yardstick(FIT_DATA, .x))
}

# time informatvie orders
orderA_params <- order_yardsticks %>%
  filter(tier == "time_informative") %>%
  transmute(Class_final = class, Order_final = order,
            ord_intercept = final_intercept, ord_slope = final_slope,
            ord_resid_sd = final_resid_sd)

# score each pair against order or class
scored_all <- yardstick_data %>%
  left_join(orderA_params,   by = c("Class_final", "Order_final")) %>%
  left_join(class_yardsticks, by = "Class_final") %>%
  mutate(
    order_scored = !is.na(ord_slope),
    reference    = if_else(order_scored, paste0("order:", Order_final), "class"),
    ref_intercept = if_else(order_scored, ord_intercept, class_intercept),
    ref_slope     = if_else(order_scored, ord_slope,     class_slope),
    ref_resid_sd  = if_else(order_scored, ord_resid_sd,  class_resid_sd)
  ) %>%
  filter(!is.na(ref_slope)) %>%  
  mutate(
    expected_log_k2p = ref_intercept + ref_slope * log_div_time,
    residual         = log_k2p - expected_log_k2p,
    std_residual     = residual / ref_resid_sd,
    p_value          = 2 * pnorm(-abs(std_residual)),
    outlier_sr_95    = abs(std_residual) > 1.96,
    outlier_sr_99    = abs(std_residual) > 2.576,
    direction        = if_else(std_residual > 0, "faster", "slower")
  ) %>%
  group_by(Class_final) %>%
  mutate(q_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# all flagged pairs 
flagged_out <- scored_all %>%
  filter(outlier_sr_95) %>%
  arrange(desc(abs(std_residual))) %>%
  dplyr::select(Class_final, Order_final, Family_final, genus_final, sp1, sp2,
                timetree_div, k2p, expected_log_k2p, residual, std_residual,
                p_value, q_fdr, reference, order_scored, direction,
                outlier_sr_95, outlier_sr_99)
write_csv(flagged_out, file.path(OUT_DIR, "flagged_outlier_pairs_by_order.csv"))

# order outlier summary
order_outlier_summary <- scored_all %>%
  filter(order_scored) %>%
  group_by(Class_final, Order_final) %>%
  summarise(n_pairs = n(),
            n_out_95 = sum(outlier_sr_95), n_out_99 = sum(outlier_sr_99),
            n_out_fdr = sum(q_fdr < FDR_ALPHA),
            pct_out_95 = 100 * mean(outlier_sr_95), .groups = "drop") %>%
  arrange(Class_final, desc(n_pairs))
write_csv(order_outlier_summary, file.path(OUT_DIR, "order_outlier_summary.csv"))

n_ord <- sum(scored_all$order_scored); n_cls <- sum(!scored_all$order_scored)
cat(sprintf(paste0("\nScored %d pairs against tier-A Order yardsticks (%d Orders) ",
                   "and %d against their Class yardstick.\n"),
            n_ord, n_distinct(orderA_params$Order_final), n_cls))
cat(sprintf("Flagged %d pairs total at |z|>1.96 (%d Order-scored, %d Class-scored).\n",
            nrow(flagged_out), sum(flagged_out$order_scored),
            sum(!flagged_out$order_scored)))

# order fit figure
ok_png <- requireNamespace("ggplot2", quietly = TRUE)
if (ok_png) {
  library(ggplot2)
  class_cols <- c(Aves = "#D55E00", Mammalia = "#56B4E9",
                  Reptilia = "#E69F00", Amphibia = "#009E73")
  scored <- yardstick_data %>% inner_join(orderA_params,
                                          by = c("Class_final", "Order_final"))
  for (cn in c("Aves", "Mammalia", "Reptilia")) {
    d <- scored %>% filter(Class_final == cn)
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(timetree_div, k2p)) +
      geom_point(alpha = 0.35, size = 0.8, colour = class_cols[[cn]], stroke = 0) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "grey20", linewidth = 0.6) +
      facet_wrap(~ Order_final, scales = "free") +
      scale_x_log10() + scale_y_log10() +
      labs(title = sprintf("Order-level genomic yardstick — %s", cn),
           x = "Divergence time (My, log)", y = "k2p (log)") +
      theme_bw(base_size = 11)
    ggsave(file.path(OUT_DIR, "figures", sprintf("divergence_by_order_%s.png", cn)),
           p, width = 11, height = 8, dpi = 350,
           device = if (capabilities("cairo")) "png" else "png")
  }
}
