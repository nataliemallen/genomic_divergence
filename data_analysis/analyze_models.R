# analyze bird and mammal modles

library(tidyverse)
library(patchwork)

RESULTS_DIR   <- "/scratch/gautschi/allen715/GD_models/results"
POSTERIOR_DIR <- file.path(RESULTS_DIR, "posteriors")
MODEL_DIR     <- file.path(RESULTS_DIR, "models")
OUT_DIR       <- file.path(RESULTS_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CLASSES     <- c("birds", "mammals")
MODEL_TYPES <- c("avg_only", "diff_only")

labels_avg <- c(
  b_z_timetree          = "Divergence Time",
  b_z_body_size_avg     = "Body Size",
  b_z_gen_time_avg      = "Generation Time",
  b_z_clutch_avg        = "Clutch Size",
  b_z_pop_density_avg   = "Population Density",
  b_z_genome_size_avg   = "Genome Size",
  b_z_gc_avg            = "GC Content",
  b_z_n50_avg           = "N50"
)

labels_diff <- c(
  b_z_timetree          = "Divergence Time",
  b_z_body_size_diff    = "Body Size",
  b_z_gen_time_diff     = "Generation Time",
  b_z_clutch_diff       = "Clutch Size",
  b_z_genome_size_diff  = "Genome Size",
  b_z_gc_diff           = "GC Content",
  b_z_n50_diff          = "N50"
)

label_map <- list(avg_only = labels_avg, diff_only = labels_diff)

# detect how many trees are present 
discover_trees <- function(class_name, model_type) {
  pattern <- sprintf("%s_tree[0-9]+_%s\\.rds", class_name, model_type)
  files   <- list.files(POSTERIOR_DIR, pattern = pattern, full.names = FALSE)
  
  tree_ids <- as.integer(
    str_extract(files, "(?<=tree)[0-9]+")
  )
  
  return(sort(tree_ids))
}

cat("=== DISCOVERING AVAILABLE POSTERIORS ===\n")
tree_inventory <- expand_grid(class = CLASSES, model_type = MODEL_TYPES) %>%
  rowwise() %>%
  mutate(
    tree_ids  = list(discover_trees(class, model_type)),
    n_trees   = length(tree_ids)
  ) %>%
  ungroup()

print(tree_inventory %>% dplyr::select(class, model_type, n_trees))

# load and pool posteriors 
summarize_posterior <- function(post, label_map_for_model) {
  fe_cols <- intersect(names(label_map_for_model), colnames(post))
  
  if (length(fe_cols) == 0) {
    warning("No matching predictor columns found — check label map keys")
    return(NULL)
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

load_and_pool <- function(class_name, model_type) {
  lmap    <- label_map[[model_type]]
  tree_ids <- discover_trees(class_name, model_type)
  
  if (length(tree_ids) == 0) {
    cat(sprintf("  No trees found for %s %s\n", class_name, model_type))
    return(NULL)
  }
  
  cat(sprintf("  Loading %s %s: %d trees\n", class_name, model_type, length(tree_ids)))
  
  # Pool all draws across trees into one dataframe
  pooled <- map_dfr(tree_ids, function(tid) {
    path <- file.path(POSTERIOR_DIR,
                      sprintf("%s_tree%04d_%s.rds", class_name, tid, model_type))
    if (!file.exists(path)) {
      warning(sprintf("Missing: %s", path))
      return(NULL)
    }
    post <- readRDS(path)
    post$tree_id <- tid
    return(post)
  })
  
  summary <- summarize_posterior(pooled, lmap)
  if (!is.null(summary)) {
    summary$class      <- class_name
    summary$model_type <- model_type
    summary$n_trees    <- length(tree_ids)
  }
  
  return(list(pooled = pooled, summary = summary, tree_ids = tree_ids))
}

cat("\n=== LOADING AND POOLING POSTERIORS ===\n")
all_posteriors <- list()

for (cn in CLASSES) {
  for (mt in MODEL_TYPES) {
    key <- paste(cn, mt, sep = "_")
    all_posteriors[[key]] <- load_and_pool(cn, mt)
  }
}

# convergence summaries
cat("\n=== CONVERGENCE SUMMARY ===\n")

convergence_all <- map_dfr(CLASSES, function(cn) {
  map_dfr(discover_trees(cn, "avg_only"), function(tid) {
    path <- file.path(MODEL_DIR,
                      sprintf("%s_tree%04d_summary.rds", cn, tid))
    if (!file.exists(path)) return(NULL)
    s <- readRDS(path)
    
    # Handle new dual-model summary structure
    avg_diag  <- s$avg_model$diagnostics
    diff_diag <- s$diff_model$diagnostics
    
    tibble(
      class           = cn,
      tree_id         = tid,
      rhat_max_avg    = avg_diag$rhat_max,
      rhat_max_diff   = diff_diag$rhat_max,
      ess_min_avg     = avg_diag$ess_min,
      ess_min_diff    = diff_diag$ess_min,
      n_divergent_avg  = avg_diag$n_divergent,
      n_divergent_diff = diff_diag$n_divergent,
      converged_avg   = avg_diag$converged,
      converged_diff  = diff_diag$converged,
      n_obs_avg       = s$avg_model$nobs,
      n_obs_diff      = s$diff_model$nobs
    )
  })
})

print(convergence_all)
write_csv(convergence_all,
          file.path(OUT_DIR, "convergence_summary_all_classes.csv"))

# credibility check for each tree - how many trees is each effect credible in?
cat("\n=== PER-TREE CREDIBILITY ===\n")

per_tree_credibility <- function(class_name, model_type) {
  lmap     <- label_map[[model_type]]
  tree_ids <- discover_trees(class_name, model_type)
  
  map_dfr(tree_ids, function(tid) {
    path <- file.path(POSTERIOR_DIR,
                      sprintf("%s_tree%04d_%s.rds", class_name, tid, model_type))
    if (!file.exists(path)) return(NULL)
    post <- readRDS(path)
    
    fe_cols <- intersect(names(lmap), colnames(post))
    
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
      mutate(tree_id    = tid,
             class      = class_name,
             model_type = model_type,
             predictor  = recode(parameter, !!!lmap))
  })
}

credibility_by_tree <- bind_rows(
  map_dfr(CLASSES, function(cn)
    map_dfr(MODEL_TYPES, function(mt)
      per_tree_credibility(cn, mt)))
)

credibility_summary <- credibility_by_tree %>%
  group_by(class, model_type, parameter, predictor) %>%
  summarize(
    n_trees_total    = n(),
    n_trees_credible = sum(credible),
    pct_credible     = 100 * mean(credible),
    median_of_medians = median(median),
    .groups = "drop"
  ) %>%
  arrange(class, model_type, desc(pct_credible))

cat("\nCredibility across trees by class and model:\n")
print(credibility_summary, n = 80)

write_csv(credibility_summary,
          file.path(OUT_DIR, "credibility_across_trees_all.csv"))

# variance components
vc_all <- map_dfr(CLASSES, function(cn) {
  map_dfr(discover_trees(cn, "avg_only"), function(tid) {
    path <- file.path(MODEL_DIR,
                      sprintf("%s_tree%04d_summary.rds", cn, tid))
    if (!file.exists(path)) return(NULL)
    s <- readRDS(path)
    
    extract_vc <- function(vc, model_label) {
      tryCatch({
        tibble(
          class      = cn,
          tree_id    = tid,
          model_type = model_label,
          sd_genus   = vc$genus$sd["Intercept", "Estimate"],
          sd_phylo   = vc$mmsp1sp2$sd["Intercept", "Estimate"],
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

vc_summary <- vc_all %>%
  group_by(class, model_type) %>%
  summarize(
    across(c(sd_genus, sd_phylo, sd_residual), 
           list(mean = mean, min = min, max = max),
           na.rm = TRUE),
    .groups = "drop"
  )

cat("\nVariance components summary:\n")
print(vc_summary)
write_csv(vc_all,     file.path(OUT_DIR, "variance_components_by_tree.csv"))
write_csv(vc_summary, file.path(OUT_DIR, "variance_components_summary.csv"))

# figures
class_colors <- c(birds = "#2980B9", mammals = "#C0392B")

# Figure 1: Forest plot — one panel per class, avg and diff side by side
plot_forest_class <- function(class_name) {
  
  avg_key  <- paste0(class_name, "_avg_only")
  diff_key <- paste0(class_name, "_diff_only")
  
  if (is.null(all_posteriors[[avg_key]]) ||
      is.null(all_posteriors[[diff_key]])) {
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
      labs(title = title, x = "Standardized effect size (beta)", 
           y = NULL, color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position  = "bottom",
            panel.grid.minor = element_blank(),
            plot.title       = element_text(face = "bold"))
  }
  
  p_avg  <- make_panel(avg_summary,  "Average model")
  p_diff <- make_panel(diff_summary, "Difference model")
  
  p_avg + p_diff +
    plot_annotation(
      title    = sprintf("%s (pooled across %d trees)",
                         str_to_title(class_name), n_trees),
      subtitle = "Points = posterior median; thick = 50% CI; thin = 95% CI"
    )
}

for (cn in CLASSES) {
  p <- plot_forest_class(cn)
  if (!is.null(p)) {
    ggsave(file.path(OUT_DIR, sprintf("%s_forest_avg_diff.pdf", cn)),
           p, width = 12, height = 6)
    ggsave(file.path(OUT_DIR, sprintf("%s_forest_avg_diff.png", cn)),
           p, width = 12, height = 6, dpi = 150)
    cat(sprintf("Saved: %s_forest_avg_diff\n", cn))
  }
}

# Figure 2: Cross-class comparison — same predictor, birds vs mammals
plot_cross_class <- function(model_type) {
  
  combined <- map_dfr(CLASSES, function(cn) {
    key <- paste(cn, model_type, sep = "_")
    if (is.null(all_posteriors[[key]])) return(NULL)
    all_posteriors[[key]]$summary %>%
      mutate(class = cn)
  })
  
  if (nrow(combined) == 0) return(NULL)
  
  predictor_order <- combined %>%
    filter(class == "birds") %>%
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
      title    = sprintf("Birds vs. Mammals: %s model",
                         gsub("_", " ", model_type)),
      subtitle = "Solid = credible; open = not credible",
      x        = "Standardized effect size (beta)",
      y        = NULL, color = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank())
}

for (mt in MODEL_TYPES) {
  p <- plot_cross_class(mt)
  if (!is.null(p)) {
    ggsave(file.path(OUT_DIR, sprintf("cross_class_%s.pdf", mt)),
           p, width = 8, height = 6)
    ggsave(file.path(OUT_DIR, sprintf("cross_class_%s.png", mt)),
           p, width = 8, height = 6, dpi = 150)
    cat(sprintf("Saved: cross_class_%s\n", mt))
  }
}

# Figure 3: Between-tree consistency per class and model
plot_tree_consistency <- function(class_name, model_type) {
  
  key     <- paste(class_name, model_type, sep = "_")
  summary <- all_posteriors[[key]]$summary
  if (is.null(summary)) return(NULL)
  
  predictor_order <- summary %>%
    arrange(desc(abs(median))) %>%
    pull(predictor) %>%
    as.character()
  
  credibility_by_tree %>%
    filter(class == class_name, model_type == model_type) %>%
    mutate(predictor = factor(as.character(predictor),
                              levels = rev(predictor_order))) %>%
    filter(!is.na(predictor)) %>%
    ggplot(aes(x = median, y = predictor, color = credible)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(position = position_jitter(height = 0.12, seed = 42),
               alpha = 0.5, size = 2) +
    geom_point(data = summary %>%
                 mutate(predictor = factor(as.character(predictor),
                                           levels = rev(predictor_order))),
               aes(x = median, y = predictor),
               color = "black", size = 4, shape = 18, inherit.aes = FALSE) +
    scale_color_manual(
      values = c("TRUE" = class_colors[[class_name]], "FALSE" = "gray60"),
      labels = c("TRUE" = "Credible in this tree", "FALSE" = "Not credible")
    ) +
    labs(
      title    = sprintf("%s %s: per-tree consistency",
                         str_to_title(class_name),
                         gsub("_", " ", model_type)),
      subtitle = "Colored dots = per-tree medians; diamond = pooled median",
      x        = "Posterior median (beta)", y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

for (cn in CLASSES) {
  for (mt in MODEL_TYPES) {
    p <- plot_tree_consistency(cn, mt)
    if (!is.null(p)) {
      ggsave(file.path(OUT_DIR, sprintf("%s_%s_tree_consistency.pdf", cn, mt)),
             p, width = 7, height = 6)
      ggsave(file.path(OUT_DIR, sprintf("%s_%s_tree_consistency.png", cn, mt)),
             p, width = 7, height = 6, dpi = 150)
    }
  }
}

# export trees
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

