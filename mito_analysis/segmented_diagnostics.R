# recompute the residual-diagnostics table (shapiro, breusch-pagan, davies breakpoint, segment slopes) from a per-pair csv

suppressWarnings(suppressMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(lmtest)  # Breusch-Pagan
  library(segmented)  # segmented regression + slope() + davies.test
}))

# config
INPUT       <- "/scratch/gautschi/allen715/GD_mito/mito_yardstick/yardstick_pairs_all.csv"
OUTPUT      <- "segmented_diagnostics_out.csv"

REAPPLY_CAP <- TRUE  # TRUE if INPUT is the uncapped frame
SEED        <- 12345
GENUS_CAP   <- 25
CLASS_LEVELS<- c("Aves", "Mammalia", "Reptilia", "Amphibia")
CLASSES     <- c("Aves", "Mammalia", "Reptilia")

# load
d <- read_csv(INPUT, show_col_types = FALSE)

# derive log columns if missing
if (!"log_k2p" %in% names(d))      d$log_k2p      <- log(d$k2p)
if (!"log_div_time" %in% names(d)) d$log_div_time <- log(d$timetree_div)

# keep valid rows; restore factor so the cap reproduces the main script's group order
d <- d %>%
  filter(is.finite(log_k2p), is.finite(log_div_time), timetree_div > 0) %>%
  mutate(Class_final = factor(Class_final, levels = CLASS_LEVELS))

# optional: reproduce the genus cap
if (REAPPLY_CAP) {
  set.seed(SEED)
  d <- d %>%
    group_by(Class_final, genus_final) %>%
    slice_sample(n = GENUS_CAP, replace = FALSE) %>%
    ungroup()
}

cat("Rows per class used for the fit (compare to n_pairs_fit):\n")
print(table(d$Class_final))

# diagnose one class
diagnose_class <- function(df_class, class_name) {
  df_class <- df_class %>% filter(Class_final == class_name)
  if (nrow(df_class) < 20) {
    cat(sprintf("  %s: only %d pairs — skipping\n", class_name, nrow(df_class)))
    return(NULL)
  }
  fit <- lm(log_k2p ~ log_div_time, data = df_class)

  sw <- tryCatch(shapiro.test(residuals(fit)), error = function(e) NULL)
  bp <- tryCatch(bptest(fit),                  error = function(e) NULL)

  seg <- tryCatch(segmented(fit, seg.Z = ~ log_div_time), error = function(e) NULL)
  seg_breakpoint_mya <- if (!is.null(seg) && !is.null(seg$psi))
    exp(seg$psi[1, "Est."]) else NA_real_
  seg_pvalue <- tryCatch(davies.test(fit, seg.Z = ~ log_div_time)$p.value,
                         error = function(e) NA_real_)

  # two joint-fit slopes: row 1 = below breakpoint, row 2 = above; col 1 = est,
  # col 2 = SE
  seg_slopes <- if (!is.null(seg))
    tryCatch(slope(seg)$log_div_time, error = function(e) NULL) else NULL
  have_two <- !is.null(seg_slopes) && is.matrix(seg_slopes) && nrow(seg_slopes) >= 2
  slope_below    <- if (have_two) seg_slopes[1, 1] else NA_real_
  slope_below_se <- if (have_two) seg_slopes[1, 2] else NA_real_
  slope_above    <- if (have_two) seg_slopes[2, 1] else NA_real_
  slope_above_se <- if (have_two) seg_slopes[2, 2] else NA_real_
  n_below <- if (!is.na(seg_breakpoint_mya))
    sum(df_class$timetree_div <  seg_breakpoint_mya) else NA_integer_
  n_above <- if (!is.na(seg_breakpoint_mya))
    sum(df_class$timetree_div >= seg_breakpoint_mya) else NA_integer_

  tibble(
    class              = class_name,
    n                  = nrow(df_class),
    shapiro_p          = if (!is.null(sw)) sw$p.value else NA_real_,
    bp_p               = if (!is.null(bp)) bp$p.value else NA_real_,
    davies_p           = seg_pvalue,
    time_dependent     = !is.na(seg_pvalue) && seg_pvalue < 0.05,
    breakpoint_mya     = seg_breakpoint_mya,
    slope_below        = slope_below,  slope_below_se = slope_below_se, n_below = n_below,
    slope_above        = slope_above,  slope_above_se = slope_above_se, n_above = n_above,
    r_squared_1line    = summary(fit)$r.squared
  )
}

out <- map_dfr(CLASSES, ~ diagnose_class(d, .x))

cat("\n=== SEGMENTED DIAGNOSTICS ===\n")
print(as.data.frame(out), digits = 3)
write_csv(out, OUTPUT)
cat(sprintf("\nWrote %s\n", OUTPUT))
cat("Note: slopes/breakpoint are only meaningful where time_dependent == TRUE;\n")
cat("      report the others as n.s. in the SI table.\n")

