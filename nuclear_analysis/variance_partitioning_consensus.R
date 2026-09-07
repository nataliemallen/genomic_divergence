### variance partitioning of log(k2p): pure-trait / shared / pure-phylogeny / residual
### (desdevises et al. 2003) from bayesian r2 of full, traits-only and phylogeny-only brms models,
### fit on the mcc consensus of the 100 pmm trees; each component carries a credible interval

library(tidyverse)
library(ape)
library(brms)
library(phytools)  # force.ultrametric, for rebuilding A from a saved tree
# also needs the phangorn package for the mcc consensus tree
# (install.packages("phangorn")); it is loaded on demand via requireNamespace()

# ── configuration ─────────────────────────────────────────────────────────────

SEED <- 12345
# point this at the directory the PMM actually wrote to (must contain
# tree_objects/ and models/<class>_scaled_df.rds)
RESULTS_DIR   <- "/scratch/gautschi/allen715/GD_models/rate_model/results"
TREE_OBJ_DIR  <- file.path(RESULTS_DIR, "tree_objects")  # the 100 per-tree pruned trees (for MCC)
SCALED_DF_DIR <- file.path(RESULTS_DIR, "models")  # <class>_scaled_df.rds (pair data)
FULLFIT_DIR   <- file.path(RESULTS_DIR, "overcorrection_diagnostics")  # cached full fits (optional reuse)
OUT_DIR       <- file.path(RESULTS_DIR, "variance_partitioning")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CLASSES     <- c("birds", "mammals")
MODEL_TYPES <- c("avg_only", "diff_only")

N_CHAINS <- 4; N_ITER <- 6000; N_WARMUP <- 2000
ADAPT_DELTA <- 0.99; MAX_TREEDEPTH <- 15
N_DRAWS_MC  <- 4000  # Monte-Carlo draws for combining the three R^2 posteriors

# each model has a different set of parameters, so priors must match:
# full      : slopes (b) + RE sd + sigma
# traits    : slopes (b) + sigma          (no random effect)
# phylo-only: RE sd + sigma               (intercept only, NO b parameter)
priors_full   <- c(prior(normal(0, 1),   class = b),
                   prior(exponential(1), class = sd),
                   prior(exponential(1), class = sigma))
priors_traits <- c(prior(normal(0, 1),   class = b),
                   prior(exponential(1), class = sigma))
priors_phylo  <- c(prior(exponential(1), class = sd),
                   prior(exponential(1), class = sigma))

get_predictors <- function(model_type) {
  if (model_type == "avg_only") {
    c("z_timetree", "z_body_size_avg", "z_gen_time_avg", "z_clutch_avg",
      "z_pop_density_avg", "z_genome_size_avg", "z_gc_avg", "z_n50_avg", "sympatric")
  } else {
    c("z_timetree", "z_body_size_diff", "z_gen_time_diff", "z_clutch_diff",
      "z_genome_size_diff", "z_gc_diff", "z_n50_diff", "sympatric")
  }
}

fit_brm <- function(form, data, A = NULL, use_re, prior) {
  args <- list(formula = as.formula(form), data = data, family = gaussian(),
               prior = prior,
               chains = N_CHAINS, cores = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
               control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
               save_pars = save_pars(all = TRUE), seed = SEED, refresh = 0)
  if (use_re) args$data2 <- list(A = A)
  do.call(brm, args)
}

r2_draws <- function(fit) as.numeric(bayes_R2(fit, summary = FALSE)[, 1])

summ <- function(x) c(median = median(x),
                      lwr = unname(quantile(x, 0.025)),
                      upr = unname(quantile(x, 0.975)))

# ── main loop ─────────────────────────────────────────────────────────────────

# representative tree = mcc consensus (maxCladeCred) of the 100 pmm trees
# the pair data (df) comes from the scaled class data, subset to the MCC tips
build_A_from_tree <- function(tree) {
  if (!is.ultrametric(tree))
    tree <- phytools::force.ultrametric(tree, method = "extend")
  vcv.phylo(tree, corr = TRUE)
}
build_mcc_tree <- function(cn) {
  tree_files <- list.files(TREE_OBJ_DIR,
                           pattern = sprintf("^%s_tree[0-9]+\\.rds$", cn),
                           full.names = TRUE)
  if (length(tree_files) == 0) return(NULL)
  trees  <- lapply(tree_files, readRDS)
  # trees can differ slightly in tips (per-tree pruning/synonyms) — reduce to the
  # set of species present in ALL trees so the clade frequencies are comparable
  common <- Reduce(intersect, lapply(trees, function(t) t$tip.label))
  trees  <- lapply(trees, function(t) drop.tip(t, setdiff(t$tip.label, common)))
  class(trees) <- "multiPhylo"
  cat(sprintf("  [%s] MCC consensus over %d trees on %d common tips\n",
              cn, length(trees), length(common)))
  if (requireNamespace("phangorn", quietly = TRUE)) {
    phangorn::maxCladeCred(trees, tree = TRUE, rooted = TRUE)
  } else {
    cat("  !! phangorn not installed (install.packages('phangorn')); ",
        "falling back to the first tree.\n", sep = "")
    trees[[1]]
  }
}
load_class_objects <- function(cn) {
  sdf_path <- file.path(SCALED_DF_DIR, sprintf("%s_scaled_df.rds", cn))
  if (!file.exists(sdf_path)) {
    cat(sprintf("!! [%s] missing scaled data %s — run the PMM first.\n", cn, sdf_path))
    return(NULL)
  }
  tree <- build_mcc_tree(cn)
  if (is.null(tree)) {
    cat(sprintf("!! [%s] no per-tree trees in %s — cannot build MCC.\n", cn, TREE_OBJ_DIR))
    return(NULL)
  }
  sdf <- readRDS(sdf_path)
  df  <- sdf %>% dplyr::filter(sp1 %in% tree$tip.label, sp2 %in% tree$tip.label)
  list(tree = tree, A = build_A_from_tree(tree), df = df)
}

partition_rows <- list()
r2_rows        <- list()

for (cn in CLASSES) {
  objs <- load_class_objects(cn)
  if (is.null(objs)) next
  tree <- objs$tree; A <- objs$A; df_full <- objs$df

  for (mt in MODEL_TYPES) {
    cat(sprintf("\n================ %s / %s ================\n", cn, mt))
    preds <- get_predictors(mt)
    df <- df_full %>%
      dplyr::select(all_of(c("log_k2p", "sp1", "sp2", "genus", preds))) %>%
      na.omit()
    cat(sprintf("  complete cases: %d\n", nrow(df)))

    pred_str <- paste(preds, collapse = " + ")
    cache <- file.path(OUT_DIR, sprintf("%s_%s_vpart_fits.rds", cn, mt))

    if (file.exists(cache)) {
      cat("  loading cached partition fits...\n")
      ff <- readRDS(cache)
      fit_full <- ff$fit_full; fit_traits <- ff$fit_traits; fit_phylo <- ff$fit_phylo
    } else {
      # reuse a cached full model if present
      full_cache <- file.path(FULLFIT_DIR, sprintf("%s_%s_fits.rds", cn, mt))
      if (file.exists(full_cache)) {
        cat("  reusing cached FULL (traits+phylo) model...\n")
        fit_full <- readRDS(full_cache)$fit_phylo
      } else {
        cat("  fitting FULL (traits + phylo) model...\n")
        fit_full <- fit_brm(sprintf("log_k2p ~ %s + (1 | mm(sp1, sp2, cov = A))", pred_str),
                            df, A, use_re = TRUE, prior = priors_full)
      }
      cat("  fitting TRAITS-ONLY model (no phylogeny)...\n")
      fit_traits <- fit_brm(sprintf("log_k2p ~ %s", pred_str), df,
                            use_re = FALSE, prior = priors_traits)
      cat("  fitting PHYLO-ONLY model (intercept + phylogeny)...\n")
      fit_phylo  <- fit_brm("log_k2p ~ 1 + (1 | mm(sp1, sp2, cov = A))", df, A,
                            use_re = TRUE, prior = priors_phylo)
      saveRDS(list(fit_full = fit_full, fit_traits = fit_traits, fit_phylo = fit_phylo), cache)
    }

    ### bayesian R^2 posteriors
    cat("  computing Bayesian R^2...\n")
    R_full   <- r2_draws(fit_full)
    R_traits <- r2_draws(fit_traits)
    R_phylo  <- r2_draws(fit_phylo)

    r2_rows[[paste(cn, mt)]] <- tibble(
      class = cn, model_type = mt,
      R2_full_med   = median(R_full),
      R2_traits_med = median(R_traits),
      R2_phylo_med  = median(R_phylo)
    )

    ### Monte-Carlo partition (three independent posteriors)
    s_full   <- sample(R_full,   N_DRAWS_MC, replace = TRUE)
    s_traits <- sample(R_traits, N_DRAWS_MC, replace = TRUE)
    s_phylo  <- sample(R_phylo,  N_DRAWS_MC, replace = TRUE)
    a <- s_full - s_phylo  # pure traits
    c <- s_full - s_traits  # pure phylogeny
    b <- s_traits + s_phylo - s_full  # shared
    d <- 1 - s_full  # residual

    for (comp in c("pure_traits", "shared", "pure_phylogeny", "residual")) {
      v <- switch(comp, pure_traits = a, shared = b, pure_phylogeny = c, residual = d)
      st <- summ(v)
      partition_rows[[paste(cn, mt, comp)]] <- tibble(
        class = cn, model_type = mt, component = comp,
        median = st["median"], lwr = st["lwr"], upr = st["upr"]
      )
    }
  }
}

if (length(partition_rows) == 0)
  stop("No class inputs were loaded, so there is nothing to partition.\n",
       "Set RESULTS_DIR to the directory the PMM wrote to. For each class the MCC\n",
       "consensus needs the per-tree trees plus the scaled data, i.e.\n",
       "  tree_objects/<class>_tree*.rds  AND  models/<class>_scaled_df.rds.")

partition <- bind_rows(partition_rows)
r2_summary <- bind_rows(r2_rows)

write_csv(partition,  file.path(OUT_DIR, "variance_partition.csv"))
write_csv(r2_summary, file.path(OUT_DIR, "component_r2_medians.csv"))

# ── FIGURE: components with 95% CrI, faceted by class/model ───────────────────

comp_order <- c("pure_phylogeny", "shared", "pure_traits", "residual")
comp_cols  <- c(pure_phylogeny = "#C0392B", shared = "#8E44AD",
                pure_traits = "#2980B9", residual = "grey70")

pf <- partition %>%
  mutate(component = factor(component, levels = comp_order),
         cell = paste(class, sub("_only", "", model_type)))

p <- ggplot(pf, aes(x = median, y = component, color = component)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lwr, xmax = upr), linewidth = 0.7, fatten = 2.5) +
  scale_color_manual(values = comp_cols) +
  facet_wrap(~ cell, ncol = 2) +
  labs(title = "Variance partitioning of genomic divergence",
       subtitle = "Desdevises et al. 2003 partition | points = median, bars = 95% CrI",
       x = "Proportion of variance", y = NULL, color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "variance_partition.png"), p, width = 10, height = 7, dpi = 150)
ggsave(file.path(OUT_DIR, "variance_partition.pdf"), p, width = 10, height = 7)

# ── print summary ─────────────────────────────────────────────────────────────

cat("\n\n=====================================================================\n")
cat("VARIANCE PARTITIONING SUMMARY (proportion of variance, median [95% CrI])\n")
cat("=====================================================================\n")
cat("\nComponent R^2 medians (the inputs):\n")
print(r2_summary, n = 50)

cat("\nPartition:\n")
partition %>%
  mutate(pretty = sprintf("%.3f [%.3f, %.3f]", median, lwr, upr)) %>%
  dplyr::select(class, model_type, component, pretty) %>%
  pivot_wider(names_from = component, values_from = pretty) %>%
  dplyr::select(class, model_type, pure_phylogeny, shared, pure_traits, residual) %>%
  print(n = 50, width = Inf)

cat("\nHOW TO READ:\n")
cat("  pure_phylogeny : divergence variation that tracks the tree but NOT the\n")
cat("                   measured traits (deep lineage effects, unmeasured biology).\n")
cat("  shared         : trait variation that is phylogenetically structured — where\n")
cat("                   body size / generation time etc. live. A large shared block\n")
cat("                   is the honest statement that these traits map onto phylogeny\n")
cat("                   and jointly predict divergence, rather than acting\n")
cat("                   independently of ancestry.\n")
cat("  pure_traits    : trait effects INDEPENDENT of phylogeny — the part a\n")
cat("                   phylogenetic model can credibly attribute to traits alone.\n")
cat("  residual       : unexplained.\n")
cat("  (shared can be slightly negative = suppression/nonlinearity; report as ~0.)\n")
cat(sprintf("\nOutputs written to: %s/\n", OUT_DIR))