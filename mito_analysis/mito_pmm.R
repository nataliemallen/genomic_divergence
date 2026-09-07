# dyadic phylogenetic mixed models for mito data
# using r/4.5.2
# libraries
library(tidyverse)
library(ape)
library(brms)
library(phytools)
library(bayesplot)
library(car)
library(readxl)
library(callr)

options(future.globals.maxSize = 16 * 1024^3)  

### parameters
RUN_CLASS <- Sys.getenv("RUN_CLASS", unset = "birds")  # birds or mammals
N_TREES      <- 100  
SEED         <- 12345
RUN_PHYLO_SIGNAL_TESTS  <- TRUE  
RUN_SENSITIVITY_ANALYSES <- TRUE  

MODEL_VARIANT <- "standard"  

N_CHAINS      <- 4  
N_ITER        <- 12000  
N_WARMUP      <- 6000  
ADAPT_DELTA   <- 0.999 
MAX_TREEDEPTH <- 15  

# bird average model needs adapt_delta 0.999 + longer warmup to converge
MAX_PAIRS_PER_GENUS <- 25

MIN_QUERY_ALIGN <- 80  
MIN_REF_ALIGN   <- 80  
MAX_DIV_TIME    <- 15  

### output directories
dir.create("results",                          showWarnings = FALSE)
dir.create("results/models",                   showWarnings = FALSE)
dir.create("results/posteriors",               showWarnings = FALSE)
dir.create("results/diagnostics",              showWarnings = FALSE)
dir.create("results/figures",                  showWarnings = FALSE)
dir.create("results/species_lists",            showWarnings = FALSE)
dir.create("results/sensitivity",              showWarnings = FALSE)
dir.create("results/checkpoints",              showWarnings = FALSE)
dir.create("results/downsampling_sensitivity", showWarnings = FALSE)
dir.create("results/tree_objects",             showWarnings = FALSE)  # per-tree pruned trees

log_file <- file(sprintf("results/analysis_log_%s.txt", RUN_CLASS), open = "wt")
sink(log_file, append = TRUE, type = "output")
sink(log_file, append = TRUE, type = "message")

# load data
setwd("/scratch/gautschi/allen715/GD_mito")
traits     <- read_excel("/scratch/gautschi/allen715/GD_models/final_dataset/Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
divergence <- read_csv("/scratch/gautschi/allen715/GD_models/final_dataset/Complete_divergence_with_mito_05-14-26.csv",
                       guess_max = Inf, show_col_types = FALSE) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(c("genome1", "genome2")), as.character))
colnames(traits)     <- make.names(colnames(traits))

# mito gc
mito_gc_tbl <- readr::read_csv(
  "/scratch/gautschi/allen715/GD_models/final_dataset/mito_gc.csv",
  show_col_types = FALSE) %>% dplyr::select(accession, mito_gc)
traits <- traits %>% dplyr::left_join(mito_gc_tbl, by = "accession")

colnames(divergence) <- make.names(colnames(divergence))

divergence <- divergence %>%
  filter(is.na(mito_k2p_flag) | mito_k2p_flag == "") %>%  # keep only QC-passing coding pairs
  mutate(
    k2p                  = suppressWarnings(as.numeric(mito_k2p)),
    mito_query_align_pct = suppressWarnings(as.numeric(mito_query_align_pct)),
    mito_ref_align_pct   = suppressWarnings(as.numeric(mito_ref_align_pct))
  )

traits <- traits %>%
  mutate(
    Genome.size           = suppressWarnings(as.numeric(Genome.size)),
    Genome.repeat.content = suppressWarnings(as.numeric(Genome.repeat.content)),
    Genome.repeat.content = na_if(Genome.repeat.content, 0),
    Contig_N50_bp         = map_dbl(Genome.Contig.N50,   parse_n50_bp),
    Scaffold_N50_bp       = map_dbl(Genome.Scaffold.N50, parse_n50_bp),
    N50_bp                = coalesce(Scaffold_N50_bp, Contig_N50_bp)
  )

# exclude flagges species
excluded_species <- traits %>%
  filter(Exclude. %in% c("yes", "Yes")) %>%
  pull(accession)
cat(sprintf("Excluding %d flagged species\n", length(excluded_species)))
divergence_clean <- divergence %>%
  filter(!genome1 %in% excluded_species,
         !genome2 %in% excluded_species)
cat(sprintf("Pairs after species exclusion: %d\n", nrow(divergence_clean)))

# filters
divergence_clean <- divergence_clean %>%
  filter(
    !is.na(mito_query_align_pct),
    !is.na(mito_ref_align_pct),
    mito_query_align_pct >= MIN_QUERY_ALIGN,
    mito_ref_align_pct   >= MIN_REF_ALIGN,
    timetree_div <= MAX_DIV_TIME | is.na(timetree_div)
  )

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

# randomize sp1 and sp2 assignment
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

dyadic_check <- divergence_clean %>%
  left_join(traits %>% dplyr::select(accession, Genome.size, Body.size..kg.),
            by = c("genome_sp1" = "accession")) %>%
  rename(gsize_sp1 = Genome.size, body_sp1 = Body.size..kg.) %>%
  left_join(traits %>% dplyr::select(accession, Genome.size, Body.size..kg.),
            by = c("genome_sp2" = "accession")) %>%
  rename(gsize_sp2 = Genome.size, body_sp2 = Body.size..kg.) %>%
  mutate(
    sp1_larger_genome        = gsize_sp1 > gsize_sp2,
    sp1_larger_body          = body_sp1 > body_sp2,
    sp1_alphabetically_first = sp1 < sp2
  )
dyadic_validation <- dyadic_check %>%
  summarize(
    pct_sp1_larger_genome = 100 * mean(sp1_larger_genome, na.rm = TRUE),
    pct_sp1_larger_body   = 100 * mean(sp1_larger_body,   na.rm = TRUE),
    pct_sp1_alpha_first   = 100 * mean(sp1_alphabetically_first, na.rm = TRUE)
  )

print(dyadic_validation)

write_csv(dyadic_validation, "results/diagnostics/dyadic_structure_validation.csv")

trait_cols <- c(
  "Body.size..kg.", "Generation.time..years.", "Clutch.size",
  "Population.Density", "area_km2", "mito_gc"
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

# pair level traits
for (col in trait_cols) {
  sp1_col <- paste0(col, "_sp1")
  sp2_col <- paste0(col, "_sp2")
  merged_data[[sp1_col]] <- suppressWarnings(as.numeric(merged_data[[sp1_col]]))
  merged_data[[sp2_col]] <- suppressWarnings(as.numeric(merged_data[[sp2_col]]))
  v1 <- merged_data[[sp1_col]]
  v2 <- merged_data[[sp2_col]]
  merged_data[[paste0(col, "_avg")]]  <- ifelse(
    is.na(v1) & is.na(v2), NA, rowMeans(cbind(v1, v2), na.rm = TRUE))
  merged_data[[paste0(col, "_diff")]] <- ifelse(
    is.na(v1) | is.na(v2), NA, abs(v1 - v2))
}

merged_data <- merged_data %>%
  mutate(
    sympatric = case_when(
      is.na(intersection_km2) ~ NA_integer_,
      intersection_km2 > 0    ~ 1L,
      TRUE                    ~ 0L
    )
  )

# filter to analysis dataset
analysis_data <- merged_data %>%
  filter(
    Class_final %in% c("Aves", "Mammalia"),
    !is.na(k2p),
    !is.na(timetree_div),
    k2p          > 0,
    timetree_div > 0
  ) %>%
  mutate(
    Class  = factor(Class_final),
    Order  = factor(Order_final),
    Family = factor(Family_final),
    genus  = factor(genus_final)
  )

print(table(analysis_data$Class_final))
# genus downsampling
downsample_genus <- function(df, max_pairs = MAX_PAIRS_PER_GENUS) {
  set.seed(SEED)
  df %>%
    group_by(genus_final) %>%
    group_modify(~ if (nrow(.x) <= max_pairs) .x else slice_sample(.x, n = max_pairs)) %>%
    ungroup()
}
genus_dist_before <- analysis_data %>%
  count(Class_final, genus_final, name = "n_pairs") %>%
  arrange(Class_final, desc(n_pairs))

print(head(genus_dist_before, 10))
analysis_data_downsampled <- downsample_genus(analysis_data)

genus_dist_after <- analysis_data_downsampled %>%
  count(Class_final, genus_final, name = "n_pairs") %>%
  arrange(Class_final, desc(n_pairs))

print(head(genus_dist_after, 10))

downsampling_report <- genus_dist_before %>%
  left_join(genus_dist_after, by = c("Class_final", "genus_final"),
            suffix = c("_before", "_after")) %>%
  mutate(n_pairs_after = replace_na(n_pairs_after, 0),
         n_removed     = n_pairs_before - n_pairs_after) %>%
  filter(n_removed > 0) %>%
  arrange(desc(n_removed))
write_csv(downsampling_report, "results/genus_downsampling_report.csv")

print(downsampling_report)
analysis_data <- analysis_data_downsampled

# by class
data_by_class <- list(
  mammals = analysis_data %>% filter(Class_final == "Mammalia"),
  birds   = analysis_data %>% filter(Class_final == "Aves")
)

for (class_name in names(data_by_class)) {
  df        <- data_by_class[[class_name]]
  n_pairs   <- nrow(df)
  n_species <- length(unique(c(df$sp1, df$sp2)))
  n_genera  <- length(unique(df$genus_final))

  n_complete <- df %>%
    filter(!is.na(Body.size..kg._avg),         !is.na(Generation.time..years._avg),
           !is.na(Clutch.size_avg),             !is.na(Population.Density_avg),
           !is.na(mito_gc_avg),                 !is.na(sympatric)) %>%
    nrow()
  cat(sprintf("%s: %d pairs, %d species, %d genera, %d complete cases\n",
              class_name, n_pairs, n_species, n_genera, n_complete))
}
data_by_class <- data_by_class[RUN_CLASS]

normalize_tree_labels <- function(tree) {
  tree$tip.label <- normalize_species(tree$tip.label)
  tree
}

# phylogenetic covariance matrix 
build_phylo_cov <- function(tree, chronos_timeout = 300) {
  if (!is.ultrametric(tree)) {
    cat("  Tree not ultrametric, running chronos...\n")
    original_tree <- tree
    tree <- tryCatch({

      callr::r(
        function(tree) ape::chronos(tree, quiet = TRUE),
        args    = list(tree = original_tree),
        timeout = chronos_timeout
      )
    }, error = function(e) {
      if (inherits(e, "callr_timeout_error") ||
          grepl("timeout|time.?out|time limit", e$message, ignore.case = TRUE)) {
        cat(sprintf("  chronos timed out after %d s — falling back to force.ultrametric\n",
                    chronos_timeout))
      } else {
        cat(sprintf("  chronos failed (%s) — falling back to force.ultrametric\n",
                    e$message))
      }
      return(NULL)
    })
    if (!is.null(tree)) {
      cat("  chronos completed successfully\n")
    }
    if (is.null(tree)) {
      tree <- tryCatch({
        phytools::force.ultrametric(original_tree, method = "extend")
      }, error = function(e) {
        cat(sprintf("  force.ultrametric also failed: %s\n", e$message))
        return(NULL)
      })
      if (!is.null(tree)) {
        cat("  Tree made ultrametric via force.ultrametric (branch extension)\n")
      }
    }
    if (is.null(tree)) {
      stop("Could not make tree ultrametric — skipping this tree")
    }
  }
  A <- vcv.phylo(tree, corr = TRUE)
  cat(sprintf("  Phylogenetic covariance matrix: %d x %d\n", nrow(A), ncol(A)))
  return(A)
}

# filter pair data to only species present in the pruned tree
prepare_brms_data <- function(df, tree) {
  keep_species <- tree$tip.label
  df_filtered  <- df %>% filter(sp1 %in% keep_species, sp2 %in% keep_species)
  n_before <- nrow(df)
  n_after  <- nrow(df_filtered)
  cat(sprintf("  Filtered data: %d -> %d pairs (%.1f%% retained)\n",
              n_before, n_after, 100 * n_after / n_before))
  return(df_filtered)
}
log and z-score transform predictors
scale_predictors <- function(df) {
  df %>%
    mutate(
      log_k2p      = log(k2p),
      log_rate     = log(k2p) - log(timetree_div),
      log_timetree = log(timetree_div),

      log_body_size_avg    = log(Body.size..kg._avg),
      log_body_size_diff   = abs(log(Body.size..kg._sp1) -
                                   log(Body.size..kg._sp2)),
      log_gen_time_avg     = log(Generation.time..years._avg),
      log_gen_time_diff    = abs(log(Generation.time..years._sp1) -
                                   log(Generation.time..years._sp2)),
      log_area_avg         = log(area_km2_avg),
      log_pop_density_avg  = log(Population.Density_avg),

      z_timetree         = scale(log_timetree)[, 1],
      z_body_size_avg    = scale(log_body_size_avg)[, 1],
      z_body_size_diff   = scale(log_body_size_diff)[, 1],
      z_gen_time_avg     = scale(log_gen_time_avg)[, 1],
      z_gen_time_diff    = scale(log_gen_time_diff)[, 1],
      z_clutch_avg       = scale(Clutch.size_avg)[, 1],
      z_clutch_diff      = scale(Clutch.size_diff)[, 1],
      z_pop_density_avg  = scale(log_pop_density_avg)[, 1],
      z_mitogc_avg       = scale(mito_gc_avg)[, 1],
      z_mitogc_diff      = scale(mito_gc_diff)[, 1]
    )
}
# synonym substitution

load_synonyms <- function() {
  bird_syn <- tryCatch(
    read_csv("/scratch/gautschi/allen715/GD_models/bird_synonyms.csv",
             show_col_types = FALSE),
    error = function(e) { warning(sprintf("Could not load bird_synonyms.csv: %s", e$message)); NULL }
  )
  mammal_syn <- tryCatch(
    read_csv("/scratch/gautschi/allen715/GD_models/mammal_synonyms.csv",
             show_col_types = FALSE),
    error = function(e) { warning(sprintf("Could not load mammal_synonyms.csv: %s", e$message)); NULL }
  )
  bind_rows(bird_syn, mammal_syn)
}
apply_synonym_substitutions <- function(tree, class_name, dataset_species,
                                        synonyms_df, tree_id, verbose = TRUE) {
  class_label <- if (class_name == "birds") "birds" else "mammals"
  syn_this_class <- synonyms_df %>%
    filter(
      class == class_label,
      !is.na(`synonym in tree`),
      str_trim(`synonym in tree`) != ""
    ) %>%
    mutate(
      dataset_name = normalize_species(`species missing in tree`),
      synonym_norm = normalize_species(`synonym in tree`)
    )
  if (nrow(syn_this_class) == 0) {
    if (verbose) cat("  No synonyms to apply for this class.\n")
    return(list(tree = tree, n_substituted = 0, substitutions = tibble()))
  }
  current_tips     <- tree$tip.label
  substitution_log <- syn_this_class %>%
    mutate(
      synonym_in_tree         = synonym_norm %in% current_tips,
      target_in_dataset       = dataset_name %in% dataset_species,
      synonym_also_in_dataset = synonym_norm %in% dataset_species,
      eligible                = synonym_in_tree & !synonym_also_in_dataset
    )
  eligible          <- substitution_log %>% filter(eligible)
  skipped_conflict  <- substitution_log %>% filter(synonym_in_tree & synonym_also_in_dataset)
  skipped_not_found <- substitution_log %>% filter(!synonym_in_tree & !is.na(synonym_norm))
  n_substituted <- 0
  if (nrow(eligible) > 0) {
    for (i in seq_len(nrow(eligible))) {
      tree$tip.label[tree$tip.label == eligible$synonym_norm[i]] <- eligible$dataset_name[i]
      n_substituted <- n_substituted + 1
    }
  }
  if (verbose) {
    cat(sprintf("  Synonym substitutions: %d applied", n_substituted))
    if (nrow(skipped_conflict)  > 0) cat(sprintf(", %d skipped (synonym also in dataset)", nrow(skipped_conflict)))
    if (nrow(skipped_not_found) > 0) cat(sprintf(", %d synonyms not found in this tree",   nrow(skipped_not_found)))
    cat("\n")
    if (n_substituted > 0) {
      cat("  Substitutions made:\n")
      for (i in seq_len(nrow(eligible)))
        cat(sprintf("    %s  →  %s\n", eligible$synonym_norm[i], eligible$dataset_name[i]))
    }
    if (nrow(skipped_conflict) > 0) {
      cat("  Skipped (synonym exists in dataset — would create duplicate):\n")
      for (i in seq_len(nrow(skipped_conflict)))
        cat(sprintf("    %s (dataset name: %s)\n",
                    skipped_conflict$synonym_norm[i], skipped_conflict$dataset_name[i]))
    }
  }
  return(list(
    tree          = tree,
    n_substituted = n_substituted,
    substitutions = eligible %>%
      dplyr::select(dataset_name, synonym_norm) %>%
      mutate(tree_id = tree_id, class = class_label),
    skipped_conflict  = skipped_conflict,
    skipped_not_found = skipped_not_found
  ))
}
summarize_synonym_substitutions <- function(log_list, out_dir = "results/species_lists") {
  dir.create(out_dir, showWarnings = FALSE)
  all_subs <- bind_rows(lapply(log_list, `[[`, "substitutions"))
  if (nrow(all_subs) == 0) {
    cat("\nNo synonym substitutions were made across any tree.\n")
    return(invisible(NULL))
  }
  sub_freq <- all_subs %>%
    count(class, dataset_name, synonym_norm, name = "n_trees_applied") %>%
    arrange(class, desc(n_trees_applied))
  write_csv(sub_freq, file.path(out_dir, "synonym_substitution_summary.csv"))
  cat("\n=== SYNONYM SUBSTITUTION SUMMARY (ALL TREES) ===\n")
  for (cl in unique(sub_freq$class)) {
    cat(sprintf("\n%s: %d unique substitutions\n", toupper(cl),
                nrow(sub_freq %>% filter(class == cl))))
    print(sub_freq %>% filter(class == cl), n = 50)
  }
  cat(sprintf("\nFull summary saved to: %s\n",
              file.path(out_dir, "synonym_substitution_summary.csv")))
  return(invisible(sub_freq))
}
# phylogenetic signal tests 
test_phylogenetic_signal <- function(df, tree, class_name) {
  cat(sprintf("\n=== PHYLOGENETIC SIGNAL: %s ===\n", class_name))
  trait_data <- df %>%
    dplyr::select(sp1, k2p) %>%
    group_by(sp1) %>%
    summarize(mean_k2p = mean(log(k2p), na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(mean_k2p))
  species_in_both <- intersect(tree$tip.label, trait_data$sp1)
  if (length(species_in_both) < 10) {
    cat("Insufficient species for phylogenetic signal test\n")
    return(NULL)
  }
  tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, species_in_both))
  trait_vec   <- trait_data %>%
    filter(sp1 %in% species_in_both) %>%
    arrange(match(sp1, tree_pruned$tip.label)) %>%
    pull(mean_k2p)
  names(trait_vec) <- tree_pruned$tip.label
  lambda_result <- tryCatch(
    phylosig(tree_pruned, trait_vec, method = "lambda", test = TRUE),
    error = function(e) { cat(sprintf("Error calculating lambda: %s\n", e$message)); NULL }
  )
  K_result <- tryCatch(
    phylosig(tree_pruned, trait_vec, method = "K", test = TRUE),
    error = function(e) { cat(sprintf("Error calculating K: %s\n", e$message)); NULL }
  )
  if (!is.null(lambda_result))
    cat(sprintf("Pagel's λ = %.3f (p = %.4f)\n", lambda_result$lambda, lambda_result$P))
  if (!is.null(K_result))
    cat(sprintf("Blomberg's K = %.3f (p = %.4f)\n", K_result$K, K_result$PIC.variance.P))
  return(list(lambda = lambda_result, K = K_result))
}
# 17) model formulas
get_model_formula <- function(class_name, model_type = "avg_only",
                               model_variant = MODEL_VARIANT) {
  if (!class_name %in% c("mammals", "birds")) {
    stop(sprintf("Unknown taxonomic group: '%s'", class_name))
  }
  if (model_variant == "rate") {
    response   <- "log_rate"
    predictors <- switch(model_type,
      avg_only  = paste("z_body_size_avg + z_gen_time_avg + z_clutch_avg +",
                        "z_pop_density_avg + z_mitogc_avg + sympatric"),
      diff_only = paste("z_body_size_diff + z_gen_time_diff + z_clutch_diff +",
                        "z_mitogc_diff + sympatric"),
      stop(sprintf("Unknown model_type: '%s'. Use 'avg_only' or 'diff_only'.", model_type))
    )
    random_effects <- "(1 | mm(sp1, sp2, cov = A))"
  } else {  # standard
    response   <- "log_k2p"
    predictors <- switch(model_type,
      avg_only  = paste("z_timetree + z_body_size_avg + z_gen_time_avg + z_clutch_avg +",
                        "z_pop_density_avg + z_mitogc_avg + sympatric"),
      diff_only = paste("z_timetree + z_body_size_diff + z_gen_time_diff + z_clutch_diff +",
                        "z_mitogc_diff + sympatric"),
      stop(sprintf("Unknown model_type: '%s'. Use 'avg_only' or 'diff_only'.", model_type))
    )
    random_effects <- "(1 | mm(sp1, sp2, cov = A))"

  }
  sprintf("%s ~ %s + %s", response, predictors, random_effects)
}

# model fitting
fit_dyadic_phylo_model <- function(df, A, class_name, model_type = "avg_only") {
  cat(sprintf("\n=== FITTING MODEL: %s (%s, variant=%s) ===\n",
              class_name, model_type, MODEL_VARIANT))
  formula_str  <- get_model_formula(class_name, model_type, MODEL_VARIANT)
  cat("\nFormula:\n  ", formula_str, "\n\n")
  response_var      <- trimws(strsplit(formula_str, "~")[[1]][1])
  predictors        <- str_extract_all(formula_str, "z_[a-z0-9_]+")[[1]]
  binary_preds      <- c("sympatric")
  binary_in_formula <- binary_preds[sapply(binary_preds, function(p) grepl(p, formula_str))]
  cols_needed <- c(response_var, "sp1", "sp2", predictors, binary_in_formula)
  if (grepl("genus", formula_str)) cols_needed <- c(cols_needed, "genus")
  df_complete <- df %>%
    dplyr::select(all_of(cols_needed)) %>%
    na.omit()
  n_complete <- nrow(df_complete)
  n_total    <- nrow(df)
  cat(sprintf("Complete cases: %d / %d (%.1f%%)\n",
              n_complete, n_total, 100 * n_complete / n_total))
  if (n_complete == 0) stop("No complete cases. Cannot fit model.")
  if (n_complete < 20) warning(sprintf("Very few complete cases (%d).", n_complete))
  cat("\nFitting model (this may take 10–30 minutes)...\n")
  fit <- brm(
    formula  = as.formula(formula_str),
    data     = df_complete,
    data2    = list(A = A),
    family   = gaussian(),
    prior    = c(
      prior(normal(0, 1),   class = b),
      prior(exponential(1), class = sd),
      prior(exponential(1), class = sigma)
    ),
    chains    = N_CHAINS,
    cores     = N_CHAINS,
    iter      = N_ITER,
    warmup    = N_WARMUP,
    thin      = 1,
    control   = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
    save_pars = save_pars(all = TRUE),
    seed      = SEED
  )
  return(fit)
}
# model diagnostics
diagnose_model <- function(fit, class_name, model_type = "avg_only") {
  cat(sprintf("\n=== DIAGNOSTICS: %s (%s) ===\n", class_name, model_type))
  cat("\nRhat values:\n")
  rhat_vals <- brms::rhat(fit)
  print(summary(rhat_vals))
  if (any(rhat_vals > 1.01, na.rm = TRUE))
    warning("Some Rhat values > 1.01, model may not have converged!")
  cat("\nDivergent transitions:\n")
  div   <- brms::nuts_params(fit, pars = "divergent__")
  n_div <- sum(subset(div, Parameter == "divergent__", Value, drop = TRUE))
  cat(sprintf("  %d divergent transitions\n", n_div))
  if (n_div > 0) warning("Divergent transitions detected!")
  cat("\nEffective sample size:\n")
  ess <- brms::neff_ratio(fit)
  print(summary(ess))
  if (any(ess < 0.1, na.rm = TRUE)) warning("Some parameters have ESS ratio < 0.1")
  return(list(
    rhat_max    = max(rhat_vals, na.rm = TRUE),
    rhat_mean   = mean(rhat_vals, na.rm = TRUE),
    ess_min     = min(ess, na.rm = TRUE),
    ess_mean    = mean(ess, na.rm = TRUE),
    n_divergent = n_div,
    converged   = max(rhat_vals, na.rm = TRUE) < 1.01
  ))
}

fit_sensitivity_models <- function(df, A, class_name) {
  if (MODEL_VARIANT != "standard") {
    cat(sprintf(
      "skipping sensitivity analysis\n",
      class_name, MODEL_VARIANT
    ))
    return(invisible(NULL))
  }
  formula_base        <- get_model_formula(class_name, "avg_only", "standard")
  formula_with_genus  <- sub("\\(1 \\| mm\\(", "(1 | genus) + (1 | mm(", formula_base)
  formula_with_family <- sub("\\(1 \\| mm\\(", "(1 | Family) + (1 | mm(", formula_base)
  predictors        <- str_extract_all(formula_base, "z_[a-z0-9_]+")[[1]]
  binary_preds      <- c("sympatric")
  binary_in_formula <- binary_preds[sapply(binary_preds, function(p) grepl(p, formula_base))]
  cols_needed <- c("log_k2p", "sp1", "sp2", "genus", "Family", predictors, binary_in_formula)
  df_complete <- df %>%
    dplyr::select(all_of(cols_needed)) %>%
    na.omit()
  fit_with_args <- function(formula) {
    brm(
      formula = as.formula(formula), data = df_complete, data2 = list(A = A),
      family  = gaussian(),
      prior   = c(prior(normal(0, 1),   class = b),
                  prior(exponential(1), class = sd),
                  prior(exponential(1), class = sigma)),
      chains  = N_CHAINS, cores = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
      control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
      seed    = SEED
    )
  }
  fit_base   <- fit_with_args(formula_base)
  fit_genus  <- fit_with_args(formula_with_genus)
  fit_family <- fit_with_args(formula_with_family)
  loo_base   <- loo(fit_base,   reloo = TRUE, reloo_args = list(chains = 1, cores = 1))
  loo_genus  <- loo(fit_genus,  reloo = TRUE, reloo_args = list(chains = 1, cores = 1))
  loo_family <- loo(fit_family, reloo = TRUE, reloo_args = list(chains = 1, cores = 1))
  comparison <- loo_compare(loo_base, loo_genus, loo_family)
  print(comparison)
  vc_base   <- VarCorr(fit_base)
  vc_genus  <- VarCorr(fit_genus)
  vc_family <- VarCorr(fit_family)

  saveRDS(
    list(fit_base   = fit_base,
         fit_genus  = fit_genus,
         fit_family = fit_family,
         loo_comparison      = comparison,
         variance_components = list(base   = vc_base,
                                    genus  = vc_genus,
                                    family = vc_family)),
    sprintf("results/sensitivity/%s_random_effects_sensitivity.rds", class_name)
  )
  return(list(fit_base   = fit_base,
              fit_genus  = fit_genus,
              fit_family = fit_family,
              comparison = comparison))
}
# load tree selection utilities 
source("tree_selection_utils.R")

# single-tree analysis 
run_single_tree_analysis <- function(tree_id, tree_config, data_by_class) {
  class_names <- names(data_by_class)
  if (tree_is_complete(tree_id, class_names)) {
    result <- list()
    for (cn in class_names) result[[cn]] <- read_checkpoint(tree_id, cn)
    return(result)
  }
  synonyms_df <- load_synonyms()
  # load trees
  loaded_trees <- load_single_trees(
    base_dir   = "/scratch/gautschi/allen715/GD_models/final_dataset/phylosets",
    mammal_n   = if ("mammals"    %in% class_names) tree_config$mammal_tree_number[tree_id]   else NULL,
    bird_file  = if ("birds"      %in% class_names) tree_config$bird_file[tree_id]            else NULL,
    bird_tree  = if ("birds"      %in% class_names) tree_config$bird_tree_index[tree_id]      else NULL,
    squam_file = if ("reptiles"   %in% class_names) tree_config$squamate_file[tree_id]        else NULL,
    squam_tree = if ("reptiles"   %in% class_names) tree_config$squamate_tree_index[tree_id]  else NULL,
    amph_file  = if ("amphibians" %in% class_names) tree_config$amphibian_file[tree_id]       else NULL,
    amph_tree  = if ("amphibians" %in% class_names) tree_config$amphibian_tree_index[tree_id] else NULL,
    verbose    = TRUE
  )
  for (name in names(loaded_trees)) {
    if (is.null(loaded_trees[[name]])) {
      cat(sprintf("  %s: NULL (PROBLEM!)\n", name))
    } else {
      cat(sprintf("  %s: OK (%d tips)\n", name, length(loaded_trees[[name]]$tip.label)))
    }
  }
  cat("\n")
  results     <- list()
  synonym_log <- list()
  for (class_name in class_names) {
    cat(sprintf("\n--- %s (tree %d) ---\n", toupper(class_name), tree_id))
    # resume
    existing <- read_checkpoint(tree_id, class_name)
    if (!is.null(existing)) {
      cat(sprintf("loading from checkpoint.\n", class_name))
      results[[class_name]] <- existing
      next
    }
    tree_key <- if (class_name == "reptiles") "squamates" else class_name
    tree     <- loaded_trees[[tree_key]]
    if (is.null(tree)) {
      warning(sprintf("No tree loaded for %s (tree_key: %s)", class_name, tree_key))
      next
    }
    cat(sprintf("  Tree loaded: %d tips\n", length(tree$tip.label)))
    tree <- normalize_tree_labels(tree)
    dataset_species <- unique(c(data_by_class[[class_name]]$sp1,
                                data_by_class[[class_name]]$sp2))
    syn_result <- apply_synonym_substitutions(
      tree            = tree,
      class_name      = class_name,
      dataset_species = dataset_species,
      synonyms_df     = synonyms_df,
      tree_id         = tree_id,
      verbose         = TRUE
    )
    tree <- syn_result$tree
    synonym_log[[length(synonym_log) + 1]] <- syn_result
    pruning_result <- prune_tree_to_dataset(
      tree       = tree,
      df         = data_by_class[[class_name]],
      class_name = class_name,
      tree_id    = tree_id,
      min_tips   = 10,
      verbose    = TRUE
    )
    if (is.null(pruning_result$tree)) next
    tree_pruned <- pruning_result$tree
    df <- prepare_brms_data(data_by_class[[class_name]], tree_pruned)
    if (nrow(df) < 10) {
      warning(sprintf("Only %d pairs after filtering for %s (tree %d) — skipping.",
                      nrow(df), class_name, tree_id))
      next
    }
    # save the exact pruned tree
    saveRDS(tree_pruned,
            sprintf("results/tree_objects/%s_tree%04d.rds", class_name, tree_id))
    A <- tryCatch(
      build_phylo_cov(tree_pruned),
      error = function(e) {
        cat(sprintf("  build_phylo_cov failed for %s tree %d: %s — skipping\n",
                    class_name, tree_id, e$message))
        return(NULL)
      }
    )
    if (is.null(A)) next
    cat(sprintf("  Ready to fit: %d pairs, %d species in tree\n",
                nrow(df), pruning_result$n_tips))
    if (tree_id == 1 && RUN_PHYLO_SIGNAL_TESTS) {
      phylo_signal <- test_phylogenetic_signal(df, tree_pruned, class_name)
      saveRDS(phylo_signal,
              sprintf("results/diagnostics/%s_phylo_signal.rds", class_name))
    }
    tryCatch({
      cat("  Fitting average-only model...\n")
      fit_avg  <- fit_dyadic_phylo_model(df, A, class_name, "avg_only")
      diag_avg <- diagnose_model(fit_avg, class_name, "avg_only")
      posterior_avg <- as_draws_df(fit_avg)
      saveRDS(posterior_avg, posterior_path(class_name, tree_id, "avg_only"))
      cat("  Fitting difference-only model...\n")
      fit_diff  <- fit_dyadic_phylo_model(df, A, class_name, "diff_only")
      diag_diff <- diagnose_model(fit_diff, class_name, "diff_only")
      posterior_diff <- as_draws_df(fit_diff)
      saveRDS(posterior_diff, posterior_path(class_name, tree_id, "diff_only"))
      model_summary <- list(
        avg_model  = list(formula = formula(fit_avg),  nobs = nobs(fit_avg),
                          diagnostics = diag_avg,
                          variance_components = VarCorr(fit_avg, summary = TRUE)),
        diff_model = list(formula = formula(fit_diff), nobs = nobs(fit_diff),
                          diagnostics = diag_diff,
                          variance_components = VarCorr(fit_diff, summary = TRUE))
      )
      saveRDS(model_summary,
              sprintf("results/models/%s_tree%04d_summary_%s.rds",
                      class_name, tree_id, MODEL_VARIANT))
      class_result <- list(
        avg_model  = list(posterior = posterior_avg,  diagnostics = diag_avg),
        diff_model = list(posterior = posterior_diff, diagnostics = diag_diff),
        synonym_log = syn_result
      )
      results[[class_name]] <- class_result
      write_checkpoint(tree_id, class_name, class_result, suffix = MODEL_VARIANT)
      cat(sprintf("  ✓ %s PMM models fitted\n", class_name))
      if (tree_id == 1 && RUN_SENSITIVITY_ANALYSES) {
        cat("  Running sensitivity analysis...\n")
        fit_sensitivity_models(df, A, class_name)
      }
      if (tree_id == 1) {
        saveRDS(list(tree = tree_pruned, A = A, df = df),
                sprintf("results/downsampling_sensitivity/%s_tree1_objects.rds", class_name))
      }
    }, error = function(e) {
      cat(sprintf("  ✗ Error fitting %s: %s\n", class_name, e$message))
      cat(sprintf("  Full error: %s\n", toString(e)))
      results[[class_name]] <<- list(error = e$message)
    })
  }
  return(c(results, list(.synonym_log = synonym_log)))
}
# multi-tree pipeline 
generate_tree_config_with_substitutes <- function(n_trees, n_reserve = 50, seed = SEED) {
  set.seed(seed)
  full_config <- generate_tree_sample_config(n_trees = n_trees + n_reserve)
  list(
    primary = full_config[1:n_trees, ],
    reserve = full_config[(n_trees + 1):(n_trees + n_reserve), ]
  )
}
run_multi_tree_analysis <- function(tree_config_list, data_by_class) {
  if (is.data.frame(tree_config_list)) {
    tree_config <- tree_config_list
    reserve     <- NULL
  } else {
    tree_config <- tree_config_list$primary
    reserve     <- tree_config_list$reserve
  }
  n_trees     <- nrow(tree_config)
  class_names <- names(data_by_class)

  n_already_done <- report_checkpoint_status(n_trees, class_names)
  if (n_already_done == n_trees) {
    cat("All trees already complete. Consolidating results.\n")
    return(consolidate_checkpoints(n_trees, class_names))
  }
  reserve_used <- 0
  all_results  <- vector("list", n_trees)
  i            <- 1
  reserve_idx  <- 1
  while (i <= n_trees) {
    if (tree_is_complete(i, class_names)) {
      cat(sprintf("  Tree %d already complete — loading from checkpoint.\n", i))
      all_results[[i]] <- consolidate_checkpoints_single(i, class_names)
      i <- i + 1
      next
    }
    cat(sprintf("\n>>> TREE %d <<<\n", i))
    result <- run_single_tree_analysis(i, tree_config, data_by_class)
    has_error <- sapply(class_names, function(cn) {
      r <- result[[cn]]
      is.null(r) ||
        isTRUE(r$skipped) ||
        (!is.null(r$avg_model) && !is.null(r$avg_model$diagnostics$error))
    })
    failed_classes <- class_names[has_error]
    if (length(failed_classes) == 0) {
      all_results[[i]] <- result
      i <- i + 1
    } else {
      cat(sprintf("  Tree %d failed for classes: %s\n",
                  i, paste(failed_classes, collapse = ", ")))
      if (!is.null(reserve) && reserve_idx <= nrow(reserve)) {
        cat(sprintf("  Substituting with reserve tree (pool index %d)\n", reserve_idx))
        tree_config[i, ] <- reserve[reserve_idx, ]
        reserve_used      <- reserve_used + 1
        reserve_idx       <- reserve_idx + 1
        sub_log_path <- "results/checkpoints/substitution_log.csv"
        sub_entry <- data.frame(
          tree_id                 = i,
          substitute_from_reserve = reserve_idx - 1,
          timestamp               = as.character(Sys.time())
        )
        if (file.exists(sub_log_path)) {
          write_csv(bind_rows(read_csv(sub_log_path, show_col_types = FALSE),
                              sub_entry), sub_log_path)
        } else {
          write_csv(sub_entry, sub_log_path)
        }
        cat(sprintf("  Retrying tree slot %d with substitute...\n", i))
      } else {
        cat(sprintf("  No reserve trees remaining — tree %d will be missing.\n", i))
        all_results[[i]] <- list(skipped = TRUE, skip_reason = "no_reserve_available")
        i <- i + 1
      }
    }
  }
  if (reserve_used > 0) {
    cat(sprintf("\n%d trees substituted from the reserve pool.\n", reserve_used))
    cat("See results/checkpoints/substitution_log.csv for details.\n")
  }
  all_synonym_logs <- unlist(
    lapply(all_results, function(tr) {
      sl <- tr[[".synonym_log"]]
      if (is.null(sl)) list() else sl
    }),
    recursive = FALSE
  )
  summarize_synonym_substitutions(all_synonym_logs)
  summarize_tree_pruning(class_names = class_names)
  return(all_results)
}
consolidate_checkpoints_single <- function(tree_id, class_names) {
  result <- list()
  for (cn in class_names) {
    cp <- read_checkpoint(tree_id, cn)
    if (!is.null(cp)) result[[cn]] <- cp
  }
  return(result)
}
# resume system
CHECKPOINT_DIR <- "results/checkpoints"
checkpoint_path <- function(tree_id, class_name, suffix = MODEL_VARIANT) {
  sprintf("%s/tree%04d_%s_%s.rds", CHECKPOINT_DIR, tree_id, class_name, suffix)
}
# path of saved posterior
posterior_path <- function(class_name, tree_id, model_type, variant = MODEL_VARIANT) {
  sprintf("results/posteriors/%s_tree%04d_%s_%s.rds",
          class_name, tree_id, model_type, variant)
}
write_checkpoint <- function(tree_id, class_name, result, suffix = MODEL_VARIANT) {
  path <- checkpoint_path(tree_id, class_name, suffix)
  saveRDS(result, path)
  cat(sprintf("  Checkpoint saved: %s\n", path))
}
read_checkpoint <- function(tree_id, class_name, suffix = MODEL_VARIANT) {
  path <- checkpoint_path(tree_id, class_name, suffix)
  if (file.exists(path)) {
    cat(sprintf("  Resuming from checkpoint: %s\n", path))
    return(readRDS(path))
  }
  return(NULL)
}
tree_is_complete <- function(tree_id, class_names, suffix = MODEL_VARIANT) {
  all(sapply(class_names, function(cn)
    file.exists(checkpoint_path(tree_id, cn, suffix))))
}
report_checkpoint_status <- function(n_trees, class_names) {
  completed <- sum(sapply(1:n_trees, tree_is_complete, class_names = class_names))
  if (completed == 0) {
    cat("  No checkpoints found — starting fresh.\n")
  } else {
    cat(sprintf("  Found checkpoints for %d / %d trees — will resume from tree %d.\n",
                completed, n_trees, completed + 1))
    partial <- sapply(1:n_trees, function(i) {
      done <- sapply(class_names, function(cn) file.exists(checkpoint_path(i, cn)))
      any(done) && !all(done)
    })
    if (any(partial))
      cat(sprintf("  Partially completed trees (will re-run missing classes): %s\n",
                  paste(which(partial), collapse = ", ")))
  }
  return(completed)
}
consolidate_checkpoints <- function(n_trees, class_names) {

  all_results <- vector("list", n_trees)
  for (i in 1:n_trees) {
    tree_result <- list()
    for (cn in class_names) {
      cp <- read_checkpoint(i, cn)
      if (!is.null(cp)) tree_result[[cn]] <- cp
    }
    all_results[[i]] <- tree_result
  }
  n_complete <- sum(sapply(all_results, function(r) length(r) == length(class_names)))
  cat(sprintf("  %d / %d trees fully complete\n", n_complete, n_trees))
  return(all_results)
}
# run models

data_by_class <- lapply(data_by_class, scale_predictors)

set.seed(SEED)
tree_config_path  <- sprintf("results/tree_config_n%d.csv",         N_TREES)
tree_reserve_path <- sprintf("results/tree_config_n%d_reserve.csv", N_TREES)
if (!file.exists(tree_config_path)) {
  config_list <- generate_tree_config_with_substitutes(N_TREES, n_reserve = 50)
  write_csv(config_list$primary, tree_config_path)
  write_csv(config_list$reserve, tree_reserve_path)
  cat(sprintf("Tree config saved: %s\n", tree_config_path))
} else {
  cat(sprintf("Loading existing tree config: %s\n", tree_config_path))
  config_list <- list(
    primary = read_csv(tree_config_path, show_col_types = FALSE),
    reserve = if (file.exists(tree_reserve_path))
                read_csv(tree_reserve_path, show_col_types = FALSE)
              else NULL
  )
}
results <- run_multi_tree_analysis(config_list, data_by_class)
save.image(file = sprintf("results/workspace_n%d_%s_%s_trees.RData",
                           N_TREES, RUN_CLASS, MODEL_VARIANT))
### calculate VIF 
compute_vif_table <- function(df, class_name, model_type = "avg_only",
                              model_variant = "standard",
                              out_dir = "results/diagnostics") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  formula_str <- get_model_formula(class_name, model_type, model_variant)
  predictors  <- str_extract_all(formula_str, "z_[a-z0-9_]+")[[1]]
  if (grepl("sympatric", formula_str)) predictors <- c(predictors, "sympatric")
  response <- if (model_variant == "rate") "log_rate" else "log_k2p"
  df_complete <- df %>%
    dplyr::select(all_of(c(response, predictors))) %>%
    na.omit()
  if (nrow(df_complete) < length(predictors) + 2) {
    warning(sprintf("Too few complete cases for VIF (%s, %s)", class_name, model_type))
    return(NULL)
  }
  ols <- lm(reformulate(predictors, response = response), data = df_complete)
  vif_vals <- car::vif(ols)
  vif_tbl <- tibble(
    class      = class_name,
    model_type = model_type,
    predictor  = names(vif_vals),
    vif        = as.numeric(vif_vals)
  ) %>% arrange(desc(vif))
  cat(sprintf("\n=== VIF: %s (%s) | n = %d ===\n",
              class_name, model_type, nrow(df_complete)))
  print(vif_tbl, n = 50)
  cat(sprintf("  Max VIF = %.2f  (>5 = moderate concern, >10 = serious)\n",
              max(vif_tbl$vif)))
  write_csv(vif_tbl,
            file.path(out_dir, sprintf("%s_%s_vif.csv", class_name, model_type)))
  return(vif_tbl)
}
cls <- RUN_CLASS
df_scaled <- data_by_class[[cls]]
vif_avg  <- compute_vif_table(df_scaled, cls, "avg_only",  MODEL_VARIANT)
vif_diff <- compute_vif_table(df_scaled, cls, "diff_only", MODEL_VARIANT)
vif_table1 <- bind_rows(vif_avg, vif_diff)
write_csv(vif_table1,
        sprintf("results/diagnostics/%s_vif_table1.csv", cls))
saveRDS(df_scaled, sprintf("results/models/%s_scaled_df.rds", cls))