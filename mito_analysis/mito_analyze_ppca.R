library(tidyverse)
library(ape)
library(brms)
library(phytools)
library(adephylo)
library(phylobase)
library(ade4)

setwd("/scratch/gautschi/allen715/GD_mito")  
TRAITS_XLSX <- "/scratch/gautschi/allen715/GD_models/final_dataset/Master_vertebrate_traits_12-08-25nma_FINAL.xlsx"

CLASSES       <- c("birds", "mammals")
MODEL_VARIANT <- "standard"
N_PPC_GLOBAL  <- 2  # number global axes
N_PPC_LOCAL   <- 1  # number local axes
PPCA_PROX     <- "Abouheif"  # tip proximity method
SEED          <- 12345
N_CHAINS <- 4; N_ITER <- 6000; N_WARMUP <- 2000
ADAPT_DELTA <- 0.99; MAX_TREEDEPTH <- 15

dir.create("results/ppca", showWarnings = FALSE, recursive = TRUE)

OVERLAP_TERMS <- "sympatric"

normalize_species <- function(x) {
  x %>% str_trim() %>% str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>% str_replace_all("[^A-Za-z ]", "") %>%
    str_to_lower() %>% str_replace_all(" ", "_")
}
parse_n50_bp <- function(x) {
  if (is.na(x)) return(NA_real_)
  x <- str_trim(x); value <- as.numeric(str_extract(x, "[0-9.]+"))
  unit <- str_extract(tolower(x), "kb|mb|gb")
  if (is.na(value) | is.na(unit)) return(NA_real_)
  value * case_when(unit == "kb" ~ 1e3, unit == "mb" ~ 1e6, unit == "gb" ~ 1e9, TRUE ~ 1)
}

# species-level trait matrix 
traits <- readxl::read_excel(TRAITS_XLSX)
colnames(traits) <- make.names(colnames(traits))

# mito gc content 
mito_gc_tbl <- readr::read_csv(
  "/scratch/gautschi/allen715/GD_models/final_dataset/mito_gc.csv",
  show_col_types = FALSE) %>% dplyr::select(accession, mito_gc)
traits <- traits %>% dplyr::left_join(mito_gc_tbl, by = "accession")

traits <- traits %>%
  mutate(
    Genome.size     = suppressWarnings(as.numeric(Genome.size)),
    Contig_N50_bp   = map_dbl(Genome.Contig.N50,   parse_n50_bp),
    Scaffold_N50_bp = map_dbl(Genome.Scaffold.N50, parse_n50_bp),
    N50_bp          = coalesce(Scaffold_N50_bp, Contig_N50_bp),
    species_full    = normalize_species(paste(genus, species))
  )

PPCA_TRAIT_SPECS <- list(
  log_body_size   = list(col = "Body.size..kg.",          fun = function(x) log(x)),
  log_gen_time    = list(col = "Generation.time..years.", fun = function(x) log(x)),
  clutch          = list(col = "Clutch.size",             fun = function(x) x),
  log_pop_density = list(col = "Population.Density",       fun = function(x) log(x)),
  mito_gc         = list(col = "mito_gc",                 fun = function(x) x)
)

# ppca functions
build_species_trait_matrix <- function(tree, traits_df) {
  tips <- tree$tip.label
  sub  <- traits_df %>% filter(species_full %in% tips)
  if (nrow(sub) == 0)
    cat("  WARNING: 0 species matched tree tips — check species_full vs tip labels\n")
  mat  <- sapply(PPCA_TRAIT_SPECS, function(s) {
    s$fun(suppressWarnings(as.numeric(sub[[s$col]])))
  })
  mat <- as.matrix(mat); rownames(mat) <- sub$species_full
  mat[!is.finite(mat)] <- NA
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  as.data.frame(mat)
}

compute_ppca <- function(tree, class_name,
                         n_global = N_PPC_GLOBAL, n_local = N_PPC_LOCAL) {
  mat <- build_species_trait_matrix(tree, traits)
  if (nrow(mat) < 10) {
    cat(sprintf("  pPCA: too few complete species (%d) for %s — skipping\n",
                nrow(mat), class_name)); return(NULL)
  }
  tree_p <- drop.tip(tree, setdiff(tree$tip.label, rownames(mat)))
  mat    <- mat[tree_p$tip.label, , drop = FALSE]

  p4d  <- phylobase::phylo4d(tree_p, tip.data = mat)
  prox <- adephylo::proxTips(tree_p, method = PPCA_PROX, normalize = "row")
  pca  <- adephylo::ppca(p4d, prox = prox, scannf = FALSE,
                         nfposi = n_global, nfnega = n_local,
                         scale = TRUE, center = TRUE)

  n_ax     <- n_global + n_local
  ax_names <- c(if (n_global > 0) paste0("ppc", seq_len(n_global)) else character(0),
                if (n_local  > 0) paste0("ppc", n_global + seq_len(n_local)) else character(0))
  ax_type  <- c(rep("global", n_global), rep("local", n_local))

  scores <- as.data.frame(pca$li)[, seq_len(n_ax), drop = FALSE]
  colnames(scores) <- ax_names; scores$species <- rownames(pca$li)
  load <- as.data.frame(pca$c1)[, seq_len(n_ax), drop = FALSE]
  colnames(load) <- ax_names;   load$trait <- rownames(pca$c1)

  moran_p <- sapply(ax_names, function(a) {
    tryCatch({
      vv <- setNames(scores[[a]], scores$species)[rownames(prox)]
      as.numeric(adephylo::abouheif.moran(vv, W = prox)$pvalue)[1]
    }, error = function(e) NA_real_)
  })

  # global-vs-local check
  cat(sprintf("\n  === pPCA AXIS CHECK: %s (%d species) ===\n", class_name, nrow(mat)))
  eig_ret <- c(if (n_global > 0) head(pca$eig, n_global) else numeric(0),
               if (n_local  > 0) tail(pca$eig, n_local)  else numeric(0))
  for (i in seq_along(ax_names)) {
    a <- ax_names[i]; typ <- ax_type[i]; p <- moran_p[[a]]
    ok <- if (is.na(p)) "? (Moran failed)" else
      if ((typ == "global" && p < 0.05) || (typ == "local" && p >= 0.05)) "OK"
      else "** UNEXPECTED — reconsider n global/local axes **"
    cat(sprintf("    %-5s [%-6s]  eig = %+.3f   Moran p = %s   %s\n",
                a, typ, eig_ret[i], ifelse(is.na(p), "NA", sprintf("%.4f", p)), ok))
  }
  cat("  (global axes want Moran p < 0.05; local axes want p >= 0.05)\n")

  list(scores = scores, loadings = load, eig = pca$eig,
       ax_type = setNames(ax_type, ax_names), moran_p = moran_p,
       n_species = nrow(mat), traits_used = names(PPCA_TRAIT_SPECS))
}

align_ppca_to_reference <- function(ppca_obj, ref_loadings) {
  if (is.null(ref_loadings)) return(ppca_obj)
  ax_names <- names(ppca_obj$ax_type)
  for (grp in c("global", "local")) {
    for (a in ax_names[ppca_obj$ax_type == grp]) {
      if (!a %in% colnames(ref_loadings)) next
      cc <- cor(ppca_obj$loadings[[a]], ref_loadings[[a]])
      if (!is.na(cc) && cc < 0) {
        ppca_obj$loadings[[a]] <- -ppca_obj$loadings[[a]]
        ppca_obj$scores[[a]]   <- -ppca_obj$scores[[a]]
      }
    }
  }
  ppca_obj
}

add_ppca_pair_predictors <- function(df, ppca_obj) {
  ax_names <- names(ppca_obj$ax_type); sc <- ppca_obj$scores
  s1 <- sc; colnames(s1) <- c(paste0(ax_names, "_s1"), "sp1")
  s2 <- sc; colnames(s2) <- c(paste0(ax_names, "_s2"), "sp2")
  df <- df %>% left_join(s1, by = "sp1") %>% left_join(s2, by = "sp2")
  for (a in ax_names) {
    df[[paste0(a, "_avg")]]  <- (df[[paste0(a, "_s1")]] + df[[paste0(a, "_s2")]]) / 2
    df[[paste0(a, "_diff")]] <- abs(df[[paste0(a, "_s1")]] - df[[paste0(a, "_s2")]])
  }
  df
}

get_ppca_predictors <- function(model_type, n_axes = N_PPC_GLOBAL + N_PPC_LOCAL) {
  suffix <- if (model_type == "ppca_avg") "_avg" else "_diff"
  ax     <- paste0("ppc", seq_len(n_axes), suffix)
  paste(c("z_timetree", ax, OVERLAP_TERMS), collapse = " + ")
}

# fit ppc axis regression
fit_ppca_regression <- function(df, class_name, model_type) {
  preds <- get_ppca_predictors(model_type)
  formula_str <- sprintf("log_k2p ~ %s", preds)
  cat(sprintf("  fitting %s: %s\n", model_type, formula_str))
  cols_needed <- c("log_k2p", "sp1", "sp2",
                   unlist(strsplit(gsub(" ", "", preds), "\\+")))
  df_complete <- df %>% dplyr::select(all_of(cols_needed)) %>% na.omit()
  if (nrow(df_complete) < 20) { warning("Too few complete cases."); return(NULL) }
  brm(formula = as.formula(formula_str), data = df_complete, family = gaussian(),
      prior = c(prior(normal(0, 1), class = b), prior(exponential(1), class = sigma)),
      chains = N_CHAINS, cores = N_CHAINS, iter = N_ITER, warmup = N_WARMUP, thin = 2,
      control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
      save_pars = save_pars(all = TRUE), seed = SEED, refresh = 0)
}

posterior_path <- function(cn, tid, mt) {
  sprintf("results/posteriors/%s_tree%04d_%s_%s.rds", cn, tid, mt, MODEL_VARIANT)
}

# test whether phylogenetic structure remains in residuals
resid_moran <- function(fit, df_used, tree) {
  r <- residuals(fit)[, "Estimate"]
  long <- bind_rows(tibble(sp = df_used$sp1, r = r),
                    tibble(sp = df_used$sp2, r = r)) %>%
    group_by(sp) %>% summarize(mr = mean(r), .groups = "drop")
  common <- intersect(tree$tip.label, long$sp)
  if (length(common) < 10) return(NA_real_)
  tr   <- drop.tip(tree, setdiff(tree$tip.label, common))
  prox <- adephylo::proxTips(tr, method = PPCA_PROX, normalize = "row")
  vv   <- setNames(long$mr[match(common, long$sp)], common)[rownames(prox)]
  tryCatch(as.numeric(adephylo::abouheif.moran(vv, W = prox)$pvalue)[1],
           error = function(e) NA_real_)
}


for (cn in CLASSES) {
  scaled_path <- sprintf("results/models/%s_scaled_df.rds", cn)
  if (!file.exists(scaled_path)) {
    cat(sprintf("\n[%s] no scaled_df (%s) — run the PMM for this class first. Skipping.\n",
                cn, scaled_path)); next
  }
  scaled_df  <- readRDS(scaled_path)
  tree_files <- list.files("results/tree_objects",
                           pattern = sprintf("^%s_tree[0-9]+\\.rds$", cn),
                           full.names = TRUE)
  tree_ids <- sort(as.integer(str_extract(basename(tree_files), "(?<=tree)[0-9]+")))

  for (tid in tree_ids) {
    tree <- readRDS(sprintf("results/tree_objects/%s_tree%04d.rds", cn, tid))
    df   <- scaled_df %>% filter(sp1 %in% tree$tip.label, sp2 %in% tree$tip.label)
    if (nrow(df) < 20) { cat(sprintf("  tree %d: too few pairs — skipping\n", tid)); next }

    done <- all(file.exists(posterior_path(cn, tid, "ppca_avg")),
                file.exists(posterior_path(cn, tid, "ppca_diff")))
    if (done) { cat(sprintf("  tree %d already complete — skipping\n", tid)); next }

    ppca_obj <- tryCatch(compute_ppca(tree, cn),
                         error = function(e) { cat(sprintf("  tree %d pPCA failed: %s\n", tid, e$message)); NULL })
    if (is.null(ppca_obj)) next

    ref_path <- sprintf("results/ppca/%s_reference_loadings.rds", cn)
    ref      <- if (file.exists(ref_path)) readRDS(ref_path) else NULL
    ppca_obj <- align_ppca_to_reference(ppca_obj, ref)
    if (is.null(ref)) saveRDS(ppca_obj$loadings, ref_path)
    saveRDS(ppca_obj, sprintf("results/ppca/%s_tree%04d_ppca.rds", cn, tid))

    df_ppca <- add_ppca_pair_predictors(df, ppca_obj)
    for (mt in c("ppca_avg", "ppca_diff")) {
      pp <- posterior_path(cn, tid, mt)
      if (file.exists(pp)) { cat(sprintf("  %s tree %d already done\n", mt, tid)); next }
      fit_p <- tryCatch(fit_ppca_regression(df_ppca, cn, mt),
                        error = function(e) { cat(sprintf("  %s failed: %s\n", mt, e$message)); NULL })
      if (is.null(fit_p)) next
      saveRDS(as_draws_df(fit_p), pp)

      df_used <- df_ppca %>%
        dplyr::select(dplyr::all_of(c("log_k2p", "sp1", "sp2",
          unlist(strsplit(gsub(" ", "", get_ppca_predictors(mt)), "\\+"))))) %>%
        na.omit()
      pm <- resid_moran(fit_p, df_used, tree)
      saveRDS(tibble(class = cn, model_type = mt, tree_id = tid, resid_moran_p = pm),
              sprintf("results/ppca/residmoran_%s_%s_tree%04d.rds", cn, mt, tid))
      cat(sprintf("    residual Moran p = %s\n",
                  ifelse(is.na(pm), "NA", sprintf("%.4f", pm))))
    }
  }
}


# residual moran summary
rm_files <- list.files("results/ppca", pattern = "^residmoran_.*\\.rds$", full.names = TRUE)
if (length(rm_files) > 0) {
  resid_moran_df <- bind_rows(lapply(rm_files, readRDS))
  write_csv(resid_moran_df, "results/ppca/ppca_residual_moran_by_tree.csv")
  summ <- resid_moran_df %>%
    group_by(class, model_type) %>%
    summarize(n_trees      = n(),
              median_p     = median(resid_moran_p, na.rm = TRUE),
              frac_p_lt_05 = mean(resid_moran_p < 0.05, na.rm = TRUE),
              .groups = "drop")
  write_csv(summ, "results/ppca/ppca_residual_moran_summary.csv")
  print(summ, n = 100)
} else {
  cat("\nNo fit found")
}

# ppca summary

N_PPC        <- N_PPC_GLOBAL + N_PPC_LOCAL
ppc_axis_tag <- c(rep("global", N_PPC_GLOBAL), rep("local", N_PPC_LOCAL))

ppca_axis_labels <- function(mt) {
  suffix  <- if (mt == "ppca_avg") "_avg" else "_diff"
  ax_keys <- paste0("b_ppc", seq_len(N_PPC), suffix)
  ax_lab  <- sprintf("pPC%d %s (%s)", seq_len(N_PPC), sub("_", "", suffix), ppc_axis_tag)
  setNames(c("Divergence Time", ax_lab, "Sympatric"),
           c("b_z_timetree", ax_keys, "b_sympatric"))
}

# axis-effect table
pool_axis_effects <- function(cn, mt) {
  files <- list.files("results/posteriors",
                      pattern = sprintf("^%s_tree[0-9]+_%s_%s\\.rds$", cn, mt, MODEL_VARIANT),
                      full.names = TRUE)
  if (length(files) == 0) return(NULL)
  lab    <- ppca_axis_labels(mt)
  pooled <- dplyr::bind_rows(lapply(files, readRDS))
  keep   <- intersect(names(lab), colnames(pooled))
  if (length(keep) == 0) return(NULL)
  pooled %>%
    dplyr::select(dplyr::all_of(keep)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "param", values_to = "v") %>%
    dplyr::group_by(param) %>%
    dplyr::summarize(median   = median(v),
                     q025     = quantile(v, 0.025),
                     q975     = quantile(v, 0.975),
                     prob_dir = pmax(mean(v > 0), mean(v < 0)),
                     credible = sign(quantile(v, 0.025)) == sign(quantile(v, 0.975)),
                     .groups  = "drop") %>%
    dplyr::mutate(class = cn, model_type = mt,
                  predictor = unname(lab[param]), n_trees = length(files))
}

axis_effects <- dplyr::bind_rows(lapply(CLASSES, function(cn)
  dplyr::bind_rows(lapply(c("ppca_avg", "ppca_diff"),
                          function(mt) pool_axis_effects(cn, mt)))))

if (nrow(axis_effects) > 0) {
  axis_effects %>%
    dplyr::select(class, model_type, predictor, median, q025, q975, credible, prob_dir) %>%
    dplyr::arrange(class, model_type, dplyr::desc(abs(median))) %>%
    print(n = 100)
  write_csv(axis_effects, "results/ppca/ppca_axis_effects.csv")
} else {
  cat("\nNo posteriors found")
}

# mean loadings
read_class_ppca <- function(cn) {
  files <- list.files("results/ppca",
                      pattern = sprintf("^%s_tree[0-9]+_ppca\\.rds$", cn),
                      full.names = TRUE)
  if (length(files) == 0) return(NULL)
  objs <- lapply(files, readRDS)
  ax   <- names(objs[[1]]$ax_type)
  load_mean <- dplyr::bind_rows(lapply(objs, function(o)
      tidyr::pivot_longer(o$loadings, dplyr::all_of(ax),
                          names_to = "axis", values_to = "loading"))) %>%
    dplyr::group_by(trait, axis) %>%
    dplyr::summarize(loading = mean(loading), .groups = "drop") %>%
    dplyr::mutate(class = cn)
  moran_sum <- dplyr::bind_rows(lapply(objs, function(o)
      tibble(axis = names(o$moran_p), moran_p = as.numeric(o$moran_p)))) %>%
    dplyr::group_by(axis) %>%
    dplyr::summarize(moran_p_median = median(moran_p, na.rm = TRUE),
                     moran_p_max    = max(moran_p,    na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(class = cn, axis_type = objs[[1]]$ax_type[axis])
  list(loadings = load_mean, moran = moran_sum)
}

meta      <- lapply(CLASSES, read_class_ppca)
loadings  <- dplyr::bind_rows(lapply(meta, function(m) if (!is.null(m)) m$loadings))
moran_tab <- dplyr::bind_rows(lapply(meta, function(m) if (!is.null(m)) m$moran))

if (nrow(loadings) > 0) {
  loadings %>% tidyr::pivot_wider(names_from = axis, values_from = loading) %>%
    dplyr::arrange(class) %>% print(n = 100)
  write_csv(loadings, "results/ppca/ppca_loadings_mean.csv")

  moran_tab <- moran_tab %>%
    dplyr::mutate(expectation = ifelse(axis_type == "global", "expect p<0.05", "expect p>=0.05"),
                  check = dplyr::case_when(
                    axis_type == "global" & moran_p_max    <  0.05 ~ "OK",
                    axis_type == "local"  & moran_p_median >= 0.05 ~ "OK",
                    TRUE ~ "** UNEXPECTED **"))
  print(moran_tab, n = 100)
  write_csv(moran_tab, "results/ppca/ppca_axis_moran.csv")

  # biplot of global axes' loadings per class
  bip <- loadings %>% dplyr::filter(axis %in% c("ppc1", "ppc2")) %>%
    tidyr::pivot_wider(names_from = axis, values_from = loading)
  if (all(c("ppc1", "ppc2") %in% names(bip))) {
    p_bip <- ggplot(bip, aes(x = ppc1, y = ppc2)) +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_vline(xintercept = 0, colour = "grey80") +
      geom_segment(aes(x = 0, y = 0, xend = ppc1, yend = ppc2),
                   arrow = arrow(length = unit(0.15, "cm")), colour = "#2980B9") +
      geom_text(aes(label = trait), size = 3, vjust = -0.4) +
      facet_wrap(~ class) +
      labs(title = "Phylogenetic PCA trait loadings (global axes)",
           subtitle = "pPC1 vs pPC2 — averaged across trees",
           x = "pPC1 loading", y = "pPC2 loading") +
      theme_minimal(base_size = 12)
    ggsave("results/ppca/ppca_biplot_loadings.png", p_bip, width = 10, height = 5, dpi = 150)
  }
}
