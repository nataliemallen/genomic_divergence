library(tidyverse)
library(patchwork)

# map_dfr removed in purrr 1.0+ — restore for compatibility
map_dfr <- function(.x, .f, ...) purrr::map(.x, .f, ...) |> bind_rows()

# ── CONFIGURATION ─────────────────────────────────────────────────────────────

RESULTS_DIR   <- "/scratch/gautschi/allen715/GD_models/rate_model/results"
POSTERIOR_DIR <- file.path(RESULTS_DIR, "posteriors")
MODEL_DIR     <- file.path(RESULTS_DIR, "models")

MODEL_VARIANT <- "standard"                    # "standard" or "rate"
CLASSES       <- c("birds", "mammals")     # both classes now complete
MODEL_TYPES   <- c("avg_only", "diff_only")

OUT_DIR <- file.path(RESULTS_DIR, "figures", MODEL_VARIANT)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── PREDICTOR LABEL MAPS ──────────────────────────────────────────────────────
# Keys must exactly match brms posterior column names (b_<predictor>).
# Values must be UNIQUE within each map — duplicate values cause a factor error
# in summarize_posterior(). sympatric is a binary predictor added to all models.

# Range overlap enters as the binary b_sympatric term (labeled "Sympatric").
labels_avg_standard <- c(
  b_z_timetree        = "Divergence Time",
  b_z_body_size_avg   = "Body Size",
  b_z_gen_time_avg    = "Generation Time",
  b_z_clutch_avg      = "Clutch Size",
  b_z_pop_density_avg = "Population Density",
  b_z_genome_size_avg = "Genome Size",
  b_z_gc_avg          = "GC Content",
  b_z_n50_avg         = "N50",
  b_sympatric         = "Sympatric"
)

labels_diff_standard <- c(
  b_z_timetree         = "Divergence Time",
  b_z_body_size_diff   = "Body Size Diff",
  b_z_gen_time_diff    = "Generation Time Diff",
  b_z_clutch_diff      = "Clutch Size Diff",
  b_z_genome_size_diff = "Genome Size Diff",
  b_z_gc_diff          = "GC Content Diff",
  b_z_n50_diff         = "N50 Diff",
  b_sympatric          = "Sympatric"
)

labels_avg_rate <- c(
  b_z_body_size_avg   = "Body Size",
  b_z_gen_time_avg    = "Generation Time",
  b_z_clutch_avg      = "Clutch Size",
  b_z_pop_density_avg = "Population Density",
  b_z_genome_size_avg = "Genome Size",
  b_z_gc_avg          = "GC Content",
  b_z_n50_avg         = "N50",
  b_sympatric         = "Sympatric"
)

labels_diff_rate <- c(
  b_z_body_size_diff   = "Body Size Diff",
  b_z_gen_time_diff    = "Generation Time Diff",
  b_z_clutch_diff      = "Clutch Size Diff",
  b_z_genome_size_diff = "Genome Size Diff",
  b_z_gc_diff          = "GC Content Diff",
  b_z_n50_diff         = "N50 Diff",
  b_sympatric          = "Sympatric"
)

# The phylogenetic PCA is fit AND summarized in the separate ppca_analysis.R;
# this script covers only the PMM trait models.
label_map <- if (MODEL_VARIANT == "rate") {
  list(avg_only = labels_avg_rate, diff_only = labels_diff_rate)
} else {
  list(avg_only = labels_avg_standard, diff_only = labels_diff_standard)
}

class_colors <- c(birds = "#2980B9", mammals = "#C0392B")

# ── FUNCTIONS — defined before any execution code ─────────────────────────────

summarize_posterior <- function(post, label_map_for_model) {
  fe_cols <- intersect(names(label_map_for_model), colnames(post))
  if (length(fe_cols) == 0) {
    warning("No matching predictor columns found — check label map keys")
    return(NULL)
  }
  # Guard against duplicate values in the label map (would crash factor())
  label_values <- unname(label_map_for_model)
  if (anyDuplicated(label_values)) {
    dupes <- label_values[duplicated(label_values)]
    stop(sprintf(
      "Duplicate values in label map: %s\nAll label map values must be unique.",
      paste(unique(dupes), collapse = ", ")
    ))
  }
  post %>%
    dplyr::select(all_of(fe_cols)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>%
    summarize(
      median        = median(value),
      mean          = mean(value),
      sd            = sd(value),
      q025          = quantile(value, 0.025),
      q975          = quantile(value, 0.975),
      q25           = quantile(value, 0.25),
      q75           = quantile(value, 0.75),
      prob_pos      = mean(value > 0),
      prob_neg      = mean(value < 0),
      credible      = sign(q025) == sign(q975),
      prob_directed = pmax(mean(value > 0), mean(value < 0)),
      direction     = ifelse(median(value) > 0, "positive", "negative"),
      .groups       = "drop"
    ) %>%
    mutate(
      predictor = recode(parameter, !!!label_map_for_model),
      predictor = factor(predictor, levels = rev(unname(label_map_for_model)))
    )
}

discover_trees <- function(class_name, model_type, variant = MODEL_VARIANT) {
  pattern <- sprintf("%s_tree[0-9]+_%s_%s\\.rds",
                     class_name, model_type, variant)
  files   <- list.files(POSTERIOR_DIR, pattern = pattern, full.names = FALSE)
  sort(as.integer(str_extract(files, "(?<=tree)[0-9]+")))
}

load_and_pool <- function(class_name, model_type, variant = MODEL_VARIANT) {
  lmap     <- label_map[[model_type]]
  tree_ids <- discover_trees(class_name, model_type, variant)
  if (length(tree_ids) == 0) {
    cat(sprintf("  No trees found for %s %s (%s)\n", class_name, model_type, variant))
    return(NULL)
  }
  cat(sprintf("  Loading %s %s [%s]: %d trees\n",
              class_name, model_type, variant, length(tree_ids)))
  pooled <- map_dfr(tree_ids, function(tid) {
    path <- file.path(POSTERIOR_DIR,
                      sprintf("%s_tree%04d_%s_%s.rds",
                              class_name, tid, model_type, variant))
    if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
    post <- readRDS(path)
    post$tree_id <- tid
    post
  })
  summary <- summarize_posterior(pooled, lmap)
  if (!is.null(summary)) {
    summary$class         <- class_name
    summary$model_type    <- model_type
    summary$model_variant <- variant
    summary$n_trees       <- length(tree_ids)
  }
  list(pooled = pooled, summary = summary, tree_ids = tree_ids)
}

per_tree_credibility <- function(class_name, model_type, variant = MODEL_VARIANT) {
  lmap     <- label_map[[model_type]]
  tree_ids <- discover_trees(class_name, model_type, variant)
  map_dfr(tree_ids, function(tid) {
    path <- file.path(POSTERIOR_DIR,
                      sprintf("%s_tree%04d_%s_%s.rds",
                              class_name, tid, model_type, variant))
    if (!file.exists(path)) return(NULL)
    post    <- readRDS(path)
    fe_cols <- intersect(names(lmap), colnames(post))
    if (length(fe_cols) == 0) return(NULL)
    post %>%
      dplyr::select(all_of(fe_cols)) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      group_by(parameter) %>%
      summarize(
        median   = median(value),
        q025     = quantile(value, 0.025),
        q975     = quantile(value, 0.975),
        credible = sign(q025) == sign(q975),
        .groups  = "drop"
      ) %>%
      mutate(
        tree_id    = tid,
        class      = class_name,
        model_type = model_type,
        predictor  = recode(parameter, !!!lmap)
      )
  })
}

plot_forest_class <- function(class_name, all_posteriors) {
  avg_key  <- paste0(class_name, "_avg_only")
  diff_key <- paste0(class_name, "_diff_only")
  if (is.null(all_posteriors[[avg_key]]) || is.null(all_posteriors[[diff_key]])) {
    cat(sprintf("Skipping forest plot for %s — missing posteriors\n", class_name))
    return(NULL)
  }
  avg_summary  <- all_posteriors[[avg_key]]$summary
  diff_summary <- all_posteriors[[diff_key]]$summary
  n_trees      <- all_posteriors[[avg_key]]$n_trees
  make_panel <- function(df, title) {
    df %>%
      arrange(desc(abs(median))) %>%
      mutate(predictor = factor(predictor, levels = rev(as.character(predictor)))) %>%
      ggplot(aes(x = median, y = predictor, color = credible)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      geom_linerange(aes(xmin = q025, xmax = q975), linewidth = 0.7) +
      geom_linerange(aes(xmin = q25,  xmax = q75),  linewidth = 1.6) +
      geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.3) +
      scale_color_manual(
        values = c("TRUE" = class_colors[[class_name]], "FALSE" = "gray50"),
        labels = c("TRUE" = "Credible", "FALSE" = "Not credible")
      ) +
      scale_x_continuous(limits = c(-0.65, 0.65)) +
      labs(title = title, x = "Standardized effect size (beta)", y = NULL, color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position  = "bottom",
            panel.grid.minor = element_blank(),
            plot.title       = element_text(face = "bold"))
  }
  p_avg  <- make_panel(avg_summary,  "Average model")
  p_diff <- make_panel(diff_summary, "Difference model")
  p_avg + p_diff +
    plot_annotation(
      title    = sprintf("%s — %s model (pooled across %d trees)",
                         str_to_title(class_name), MODEL_VARIANT, n_trees),
      subtitle = "Points = posterior median; thick = 50% CI; thin = 95% CI"
    )
}

plot_cross_class <- function(model_type, all_posteriors) {
  combined <- map_dfr(CLASSES, function(cn) {
    key <- paste(cn, model_type, sep = "_")
    if (is.null(all_posteriors[[key]])) return(NULL)
    all_posteriors[[key]]$summary %>% mutate(class = cn)
  })
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  predictor_order <- combined %>%
    filter(class == CLASSES[1]) %>%
    arrange(desc(abs(median))) %>%
    pull(predictor) %>%
    as.character()
  combined %>%
    mutate(
      predictor = factor(as.character(predictor), levels = rev(predictor_order)),
      class     = str_to_title(class)
    ) %>%
    filter(!is.na(predictor)) %>%
    ggplot(aes(x = median, y = predictor, color = class, shape = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_linerange(aes(xmin = q025, xmax = q975),
                   position = position_dodge(0.5), linewidth = 0.6) +
    geom_linerange(aes(xmin = q25, xmax = q75),
                   position = position_dodge(0.5), linewidth = 1.4) +
    geom_point(position = position_dodge(0.5), size = 3) +
    scale_color_manual(values = c(Birds = "#2980B9", Mammals = "#C0392B")) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Credible", "FALSE" = "Not credible")) +
    scale_x_continuous(limits = c(-0.8, 0.8)) +
    labs(
      title    = sprintf("Birds vs. Mammals: %s model (%s)",
                         gsub("_", " ", model_type), MODEL_VARIANT),
      subtitle = "Solid = credible; open = not credible",
      x = "Standardized effect size (beta)", y = NULL, color = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

plot_tree_consistency <- function(class_name, model_type_arg,
                                   credibility_by_tree, all_posteriors) {
  key     <- paste(class_name, model_type_arg, sep = "_")
  summary <- all_posteriors[[key]]$summary
  if (is.null(summary)) return(NULL)
  predictor_order <- summary %>%
    arrange(desc(abs(median))) %>%
    pull(predictor) %>%
    as.character()
  credibility_by_tree %>%
    filter(class == class_name, model_type == .env$model_type_arg) %>%
    mutate(predictor = factor(as.character(predictor),
                              levels = rev(predictor_order))) %>%
    filter(!is.na(predictor)) %>%
    ggplot(aes(x = median, y = predictor, color = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(position = position_jitter(height = 0.12, seed = 42),
               alpha = 0.5, size = 2) +
    geom_point(
      data = summary %>%
        mutate(predictor = factor(as.character(predictor),
                                   levels = rev(predictor_order))),
      aes(x = median, y = predictor),
      color = "black", size = 4, shape = 18, inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c("TRUE" = class_colors[[class_name]], "FALSE" = "gray60"),
      labels = c("TRUE" = "Credible in this tree", "FALSE" = "Not credible")
    ) +
    labs(
      title    = sprintf("%s %s: per-tree consistency (%s)",
                         str_to_title(class_name),
                         gsub("_", " ", model_type_arg),
                         MODEL_VARIANT),
      subtitle = "Colored dots = per-tree medians; diamond = pooled median",
      x = "Posterior median (beta)", y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# ── EXECUTION ─────────────────────────────────────────────────────────────────

cat("=== DISCOVERING AVAILABLE POSTERIORS ===\n")
tree_inventory <- expand_grid(class = CLASSES, model_type = MODEL_TYPES) %>%
  rowwise() %>%
  mutate(
    tree_ids = list(discover_trees(class, model_type)),
    n_trees  = length(tree_ids)
  ) %>%
  ungroup()
print(tree_inventory %>% dplyr::select(class, model_type, n_trees))

cat("\n=== LOADING AND POOLING POSTERIORS ===\n")
all_posteriors <- list()
for (cn in CLASSES) {
  for (mt in MODEL_TYPES) {
    key <- paste(cn, mt, sep = "_")
    all_posteriors[[key]] <- load_and_pool(cn, mt)
  }
}

for (key in names(all_posteriors)) {
  if (!is.null(all_posteriors[[key]])) {
    cat(sprintf("\n=== %s ===\n", toupper(key)))
    all_posteriors[[key]]$summary %>%
      dplyr::select(predictor, median, q025, q975, prob_directed, direction, credible) %>%
      arrange(desc(abs(median))) %>%
      print()
  }
}

cat("\n=== CONVERGENCE SUMMARY ===\n")
convergence_all <- map_dfr(CLASSES, function(cn) {
  map_dfr(discover_trees(cn, "avg_only"), function(tid) {
    path <- file.path(MODEL_DIR,
                      sprintf("%s_tree%04d_summary_%s.rds", cn, tid, MODEL_VARIANT))
    if (!file.exists(path)) {
      path <- file.path(MODEL_DIR, sprintf("%s_tree%04d_summary.rds", cn, tid))
    }
    if (!file.exists(path)) return(NULL)
    s <- readRDS(path)
    avg_diag  <- s$avg_model$diagnostics
    diff_diag <- s$diff_model$diagnostics
    tibble(
      class            = cn,
      tree_id          = tid,
      rhat_max_avg     = avg_diag$rhat_max,
      rhat_max_diff    = diff_diag$rhat_max,
      ess_min_avg      = avg_diag$ess_min,
      ess_min_diff     = diff_diag$ess_min,
      n_divergent_avg  = avg_diag$n_divergent,
      n_divergent_diff = diff_diag$n_divergent,
      converged_avg    = avg_diag$converged,
      converged_diff   = diff_diag$converged,
      n_obs_avg        = s$avg_model$nobs,
      n_obs_diff       = s$diff_model$nobs
    )
  })
})
print(convergence_all)
write_csv(convergence_all,
          file.path(OUT_DIR, "convergence_summary_all_classes.csv"))

cat("\n=== PER-TREE CREDIBILITY ===\n")
credibility_by_tree <- bind_rows(
  map_dfr(CLASSES,     function(cn)
    map_dfr(MODEL_TYPES, function(mt)
      per_tree_credibility(cn, mt)))
)

if (nrow(credibility_by_tree) > 0) {
  credibility_summary <- credibility_by_tree %>%
    group_by(class, model_type, parameter, predictor) %>%
    summarize(
      n_trees_total     = n(),
      n_trees_credible  = sum(credible),
      pct_credible      = 100 * mean(credible),
      median_of_medians = median(median),
      .groups = "drop"
    ) %>%
    arrange(class, model_type, desc(pct_credible))
  cat("\nCredibility across trees:\n")
  print(credibility_summary, n = 80)
  write_csv(credibility_summary,
            file.path(OUT_DIR, "credibility_across_trees_all.csv"))
} else {
  cat("No credibility data found — check posterior file paths\n")
  credibility_summary <- tibble()
}

cat("\n=== VARIANCE COMPONENTS ACROSS TREES ===\n")
vc_all <- map_dfr(CLASSES, function(cn) {
  map_dfr(discover_trees(cn, "avg_only"), function(tid) {
    path <- file.path(MODEL_DIR,
                      sprintf("%s_tree%04d_summary_%s.rds", cn, tid, MODEL_VARIANT))
    if (!file.exists(path)) {
      path <- file.path(MODEL_DIR, sprintf("%s_tree%04d_summary.rds", cn, tid))
    }
    if (!file.exists(path)) return(NULL)
    s <- readRDS(path)
    extract_vc <- function(vc, model_label) {
      # genus is no longer a random effect in either model variant —
      # phylogeny (mmsp1sp2) captures genus-level structure directly.
      tryCatch({
        tibble(
          class       = cn,
          tree_id     = tid,
          model_type  = model_label,
          sd_phylo    = vc$mmsp1sp2$sd["Intercept", "Estimate"],
          sd_residual = as.numeric(vc$residual__$sd[1, "Estimate"])
        )
      }, error = function(e) NULL)
    }
    bind_rows(
      extract_vc(s$avg_model$variance_components,  "avg_only"),
      extract_vc(s$diff_model$variance_components, "diff_only")
    )
  })
})

if (!is.null(vc_all) && nrow(vc_all) > 0) {
  vc_summary <- vc_all %>%
    group_by(class, model_type) %>%
    summarize(
      across(c(sd_phylo, sd_residual),
             list(mean = ~ mean(.x, na.rm = TRUE),
                  min  = ~ min(.x,  na.rm = TRUE),
                  max  = ~ max(.x,  na.rm = TRUE))),
      .groups = "drop"
    )
  cat("\nVariance components summary:\n")
  print(vc_summary)
  write_csv(vc_all,     file.path(OUT_DIR, "variance_components_by_tree.csv"))
  write_csv(vc_summary, file.path(OUT_DIR, "variance_components_summary.csv"))
} else {
  cat("No variance component data found\n")
}

cat("\n=== GENERATING FIGURES ===\n")

for (cn in CLASSES) {
  p <- plot_forest_class(cn, all_posteriors)
  if (!is.null(p)) {
    ggsave(file.path(OUT_DIR, sprintf("%s_forest_avg_diff.pdf", cn)),
           p, width = 12, height = 6)
    ggsave(file.path(OUT_DIR, sprintf("%s_forest_avg_diff.png", cn)),
           p, width = 12, height = 6, dpi = 150)
    cat(sprintf("Saved: %s_forest_avg_diff\n", cn))
  }
}

if (length(CLASSES) > 1) {
  for (mt in MODEL_TYPES) {
    p <- plot_cross_class(mt, all_posteriors)
    if (!is.null(p)) {
      ggsave(file.path(OUT_DIR, sprintf("cross_class_%s.pdf", mt)),
             p, width = 8, height = 6)
      ggsave(file.path(OUT_DIR, sprintf("cross_class_%s.png", mt)),
             p, width = 8, height = 6, dpi = 150)
      cat(sprintf("Saved: cross_class_%s\n", mt))
    }
  }
}

if (nrow(credibility_by_tree) > 0) {
  for (cn in CLASSES) {
    for (mt in MODEL_TYPES) {
      p <- plot_tree_consistency(cn, mt, credibility_by_tree, all_posteriors)
      if (!is.null(p)) {
        ggsave(file.path(OUT_DIR, sprintf("%s_%s_tree_consistency.pdf", cn, mt)),
               p, width = 7, height = 6)
        ggsave(file.path(OUT_DIR, sprintf("%s_%s_tree_consistency.png", cn, mt)),
               p, width = 7, height = 6, dpi = 150)
        cat(sprintf("Saved: %s_%s_tree_consistency\n", cn, mt))
      }
    }
  }
}

# ── STEP 6b: RAW PREDICTOR SCATTERPLOTS ───────────────────────────────────────
# Moved here from the model scripts so plots can be regenerated WITHOUT refitting.
# Each model script saves its scaled pair-level data frame via
#   saveRDS(df_scaled, "results/models/<class>_scaled_df.rds")
# Predictors are derived from label_map (keys are "b_z_*"; strip "b_" -> df cols),
# so this needs no get_model_formula() (which lives only in the model scripts).
 
cat("\n=== RAW PREDICTOR SCATTERPLOTS ===\n")
 
plot_raw_relationships <- function(df, class_name, model_type,
                                   out_dir = OUT_DIR,
                                   response = "k2p") {   # raw k2p; log10 axis applied below
  predictors <- sub("^b_", "", names(label_map[[model_type]]))
  predictors <- intersect(predictors, colnames(df))
  if (length(predictors) == 0 || !response %in% colnames(df)) {
    cat(sprintf("  %s %s: no matching columns in scaled df — skipping\n",
                class_name, model_type))
    return(NULL)
  }
 
  pretty <- c(
    z_timetree = "Divergence time", z_body_size_avg = "Body size (avg)",
    z_gen_time_avg = "Generation time (avg)", z_clutch_avg = "Clutch size (avg)",
    z_pop_density_avg = "Population density (avg)", z_genome_size_avg = "Genome size (avg)",
    z_gc_avg = "GC content (avg)", z_n50_avg = "N50 (avg)",
    z_body_size_diff = "Body size (diff)", z_gen_time_diff = "Generation time (diff)",
    z_clutch_diff = "Clutch size (diff)", z_genome_size_diff = "Genome size (diff)",
    z_gc_diff = "GC content (diff)", z_n50_diff = "N50 (diff)"
  )
 
  long <- df %>%
    dplyr::select(all_of(c(response, predictors))) %>%
    pivot_longer(all_of(predictors), names_to = "predictor", values_to = "x") %>%
    filter(!is.na(x), !is.na(.data[[response]])) %>%
    mutate(predictor = recode(predictor, !!!pretty))
 
  p <- ggplot(long, aes(x = x, y = .data[[response]])) +
    geom_point(alpha = 0.3, size = 0.9, color = "gray30") +
    geom_smooth(method = "lm", se = TRUE, color = "#C0392B", linewidth = 0.9) +
    facet_wrap(~ predictor, scales = "free_x") +
    # k2p is positive but tiny (~0.002-0.15) and spans ~100x, so a linear axis
    # collapses every point to the bottom. A log axis spreads them out; the
    # explicit breaks/labels keep the tick values positive and readable.
    # (ggplot applies the scale BEFORE fitting, so the red lm line is a true
    #  best fit on this axis.) For a plain linear axis instead, delete this line.
    scale_y_log10(breaks = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1),
                  labels = c("0.002", "0.005", "0.01", "0.02", "0.05", "0.1")) +
    labs(
      title    = sprintf("%s: k2p divergence vs. predictors (%s model)",
                         str_to_title(class_name), model_type),
      subtitle = "Red = linear best fit (with 95% band). y = k2p on a log axis (positive values; spans ~100x).",
      x = "Standardized predictor (z-score)",
      y = "k2p divergence"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
 
  ggsave(file.path(out_dir, sprintf("%s_%s_raw_scatter.png", class_name, model_type)),
         p, width = 11, height = 8, dpi = 150)
  ggsave(file.path(out_dir, sprintf("%s_%s_raw_scatter.pdf", class_name, model_type)),
         p, width = 11, height = 8)
  cat(sprintf("  Saved: %s_%s_raw_scatter\n", class_name, model_type))
  return(p)
}
 
# Locate each class's scaled df. The model scripts setwd() into .../rate_model,
# so the file may live under a different results/ tree than MODEL_DIR. Try a few
# candidate locations and use the first that exists.
find_scaled_df <- function(cn) {
  cands <- c(
    file.path(MODEL_DIR, sprintf("%s_scaled_df.rds", cn)),
    sprintf("/scratch/gautschi/allen715/GD_models/rate_model/results/models/%s_scaled_df.rds", cn),
    sprintf("/scratch/gautschi/allen715/GD_models/standard/results/models/%s_scaled_df.rds", cn)
  )
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) NA_character_ else hit[1]
}
 
for (cn in CLASSES) {
  dfp <- find_scaled_df(cn)
  if (is.na(dfp)) {
    cat(sprintf("  No scaled df found for %s — run its model script once (it saves\n", cn))
    cat(sprintf("    results/models/%s_scaled_df.rds before figures) — skipping plots.\n", cn))
    next
  }
  cat(sprintf("  %s: reading %s\n", cn, dfp))
  df_scaled <- readRDS(dfp)
  for (mt in MODEL_TYPES) plot_raw_relationships(df_scaled, cn, mt)
}
 
# ── STEP 7: EXPORT ALL TABLES ─────────────────────────────────────────────────
 
all_summaries <- bind_rows(
  map(all_posteriors, ~ if (!is.null(.x)) .x$summary else NULL)
)
 
write_csv(all_summaries,
          file.path(OUT_DIR, "pooled_effects_all_classes_models.csv"))
write_csv(credibility_summary,
          file.path(OUT_DIR, "credibility_summary_all_classes_models.csv"))
 
cat("\n=== COMPLETE ===\n")
cat(sprintf("Figures saved to: %s\n", OUT_DIR))
cat(sprintf("Tables saved to: %s\n", OUT_DIR))
