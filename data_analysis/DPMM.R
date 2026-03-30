# dyadic phylogenetic mixed model pipeline for birds or mammals

library(ape)
library(brms)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(readxl)
library(readr)
library(phytools)  
library(tidyverse)
library(bayesplot)
library(future)
library(furrr)
library(car)  

RUN_CLASS <- "birds" 

# analysis parameters
N_TREES <- 10  
SEED <- 12345
USE_PARALLEL <- FALSE  
N_CORES <- 64  

# genus downsampling parameters
MAX_PAIRS_PER_GENUS <- 25 
DOWNSAMPLE_GENERA <- c("Falco")  # always downsample

# options
RUN_PHYLO_SIGNAL_TESTS <- TRUE  # calculate Pagel's λ and Blomberg's K
RUN_MISSING_DATA_ANALYSIS <- TRUE  # analyze alignment failure patterns
RUN_SENSITIVITY_ANALYSES <- TRUE  # compare genus vs family vs phylogeny-only models
SEPARATE_AVG_DIFF_MODELS <- TRUE  # fit separate models for average vs difference effects

# MCMC parameters (Hamiltonian Monte Carlo sampling)
N_CHAINS <- 4  # number of independent chains
N_ITER <- 6000  # total iterations per chain
N_WARMUP <- 2000  # warmup/burn-in iterations
ADAPT_DELTA <- 0.99  # target acceptance rate (higher = more accurate, slower)
MAX_TREEDEPTH <- 15  # max tree depth for sampler

# genome alignment and divergence time filters
MIN_QUERY_ALIGN <- 80   # minimum query alignment % 
MIN_REF_ALIGN   <- 80   # minimum reference alignment %
MAX_DIV_TIME    <- 15   # maximum divergence time in My

# output directories
dir.create("results", showWarnings = FALSE)
dir.create("results/models", showWarnings = FALSE)
dir.create("results/posteriors", showWarnings = FALSE)
dir.create("results/diagnostics", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/species_lists", showWarnings = FALSE)
dir.create("results/missing_data", showWarnings = FALSE)
dir.create("results/sensitivity", showWarnings = FALSE)

# logging setup
log_file <- file(sprintf("results/analysis_log_%s.txt", RUN_CLASS), open = "wt")
sink(log_file, append = TRUE, type = "output")  
sink(log_file, append = TRUE, type = "message")

cat("=================================================\n")
cat("GENOMIC DIVERGENCE ANALYSIS\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat(sprintf("N_TREES: %d\n", N_TREES))
cat(sprintf("USE_PARALLEL: %s\n", USE_PARALLEL))
cat("=================================================\n\n")


# 1) load data

setwd("/scratch/gautschi/allen715/GD_models")
traits <- read_excel("/scratch/gautschi/allen715/GD_models/final_dataset/Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
divergence <- read_csv("/scratch/gautschi/allen715/GD_models/final_dataset/Complete_divergence_03-09-26.csv")

colnames(traits) <- make.names(colnames(traits))
colnames(divergence) <- make.names(colnames(divergence))

cat(sprintf("Loaded %d species with traits\n", nrow(traits)))
cat(sprintf("Loaded %d species pairs with divergence data\n", nrow(divergence)))

# 2) parse genome metrics

parse_n50_bp <- function(x) {
  if (is.na(x)) return(NA_real_)
  x <- str_trim(x)
  value <- as.numeric(str_extract(x, "[0-9.]+"))
  unit  <- str_extract(tolower(x), "kb|mb|gb")
  if (is.na(value) | is.na(unit)) return(NA_real_)
  value * case_when(
    unit == "kb" ~ 1e3,
    unit == "mb" ~ 1e6,
    unit == "gb" ~ 1e9,
    TRUE ~ 1
  )
}

traits <- traits %>%
  mutate(
    # convert genome size to numeric 
    Genome.size = suppressWarnings(as.numeric(Genome.size)),
    # convert repeat content to numeric, treat 0 as missing
    Genome.repeat.content = suppressWarnings(as.numeric(Genome.repeat.content)),
    Genome.repeat.content = na_if(Genome.repeat.content, 0),
    # parse N50 values (convert "3.2 Mb" → 3200000)
    Contig_N50_bp   = map_dbl(Genome.Contig.N50, parse_n50_bp),
    Scaffold_N50_bp = map_dbl(Genome.Scaffold.N50, parse_n50_bp),
    # use scaffold N50 if available, otherwise contig N50
    N50_bp = coalesce(Scaffold_N50_bp, Contig_N50_bp)
  )

# 3) exclude flagged species

excluded_species <- traits %>%
  filter(Exclude. %in% c("yes", "Yes")) %>%
  pull(accession)

cat(sprintf("Excluding %d flagged species\n", length(excluded_species)))

divergence_clean <- divergence %>%
  filter(!genome1 %in% excluded_species,
         !genome2 %in% excluded_species)

cat(sprintf("Pairs after species exclusion: %d\n", nrow(divergence_clean)))


# 4) alignment and divergence time filters

divergence_clean <- divergence_clean %>%
  filter(
    # must have alignment data (excludes pipeline failures)
    !is.na(query_alignment_percent),
    !is.na(ref_alignment_percent),
    # both genomes must align sufficiently
    query_alignment_percent >= MIN_QUERY_ALIGN,
    ref_alignment_percent   >= MIN_REF_ALIGN,
    # divergence time cap for birds/mammals comparability
    timetree_div <= MAX_DIV_TIME | is.na(timetree_div)
  )

cat(sprintf("Pairs after alignment filters (>=%d%% both genomes): %d\n",
            MIN_QUERY_ALIGN, nrow(divergence_clean)))

# report what was removed and why
cat(sprintf("  Removed due to missing alignment data: %d\n",
            sum(is.na(divergence$query_alignment_percent) & !divergence$genome1 %in% excluded_species)))
cat(sprintf("  Removed due to low alignment (<80%%): %d\n",
            sum(!is.na(divergence$query_alignment_percent) &
                  (divergence$query_alignment_percent < MIN_QUERY_ALIGN |
                     divergence$ref_alignment_percent < MIN_REF_ALIGN),
                na.rm = TRUE)))
cat(sprintf("  Removed due to divergence time >%d My: %d\n",
            MAX_DIV_TIME,
            sum(!is.na(divergence$timetree_div) & divergence$timetree_div > MAX_DIV_TIME)))


# 5) name normaliation

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
         genus_normalized = str_to_lower(genus))

divergence_clean <- divergence_clean %>%
  mutate(
    species1_norm = normalize_species(species1),
    species2_norm = normalize_species(species2)
  )

# 6) randomize sp1 and sp2 assignment

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

# 7) check that species randomization worked

dyadic_check <- divergence_clean %>%
  left_join(traits %>% dplyr::select(accession, Genome.size, Body.size..kg.),
            by = c("genome_sp1" = "accession")) %>%
  rename(gsize_sp1 = Genome.size, body_sp1 = Body.size..kg.) %>%
  left_join(traits %>% dplyr::select(accession, Genome.size, Body.size..kg.),
            by = c("genome_sp2" = "accession")) %>%
  rename(gsize_sp2 = Genome.size, body_sp2 = Body.size..kg.) %>%
  mutate(
    sp1_larger_genome    = gsize_sp1 > gsize_sp2,
    sp1_larger_body      = body_sp1 > body_sp2,
    sp1_alphabetically_first = sp1 < sp2
  )

dyadic_validation <- dyadic_check %>%
  summarize(
    pct_sp1_larger_genome = 100 * mean(sp1_larger_genome, na.rm = TRUE),
    pct_sp1_larger_body   = 100 * mean(sp1_larger_body,   na.rm = TRUE),
    pct_sp1_alpha_first   = 100 * mean(sp1_alphabetically_first, na.rm = TRUE)
  )

cat("Dyadic structure validation (should all be ~50%):\n")
print(dyadic_validation)

if (any(abs(dyadic_validation - 50) > 10, na.rm = TRUE)) {
  warning("Dyadic structure may not be random. Check sp1/sp2 assignment.")
}

write_csv(dyadic_validation, "results/diagnostics/dyadic_structure_validation.csv")

# 8) merge traits to pairs

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
    # genus_final will be lowercase because normalize_species() lowercased it
    # use the raw genus from traits (not normalized) to preserve case for readability,
    # but coerce to lowercase for consistency with sp1/sp2 name matching
    genus_final  = str_to_lower(coalesce(genus_sp1, genus_sp2))
  )


# 9) calculate pair-level traits

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


# 10) filter to analysis dataset

analysis_data <- merged_data %>%
  filter(
    Class_final %in% c("Aves", "Mammalia"),
    !is.na(k2p),
    !is.na(timetree_div),
    k2p > 0,
    timetree_div > 0
  ) %>%
  mutate(
    Class  = factor(Class_final),
    Order  = factor(Order_final),
    Family = factor(Family_final),
    genus  = factor(genus_final)
  )

cat(sprintf("\nFinal analysis pairs: %d\n", nrow(analysis_data)))
cat("By taxonomic group:\n")
print(table(analysis_data$Class_final))

# 11) genus downsampling

downsample_genus <- function(df, max_pairs = MAX_PAIRS_PER_GENUS) {
  set.seed(SEED)
  df %>%
    group_by(genus_final) %>%
    group_modify(~ {
      if (nrow(.x) <= max_pairs) .x else slice_sample(.x, n = max_pairs)
    }) %>%
    ungroup()
}

genus_dist_before <- analysis_data %>%
  count(Class_final, genus_final, name = "n_pairs") %>%
  arrange(Class_final, desc(n_pairs))

cat("\nTop 10 genera by pair count (before downsampling):\n")
print(head(genus_dist_before, 10))

analysis_data_downsampled <- downsample_genus(analysis_data)

# Validate downsampling 
cat("\n=== DOWNSAMPLING VALIDATION ===\n")
cat(sprintf("Falco pairs before: %d\n",
            sum(analysis_data$genus_final == "falco" &
                  analysis_data$Class_final == "Aves")))
cat(sprintf("Falco pairs after:  %d (should be <= %d)\n",
            sum(analysis_data_downsampled$genus_final == "falco" &
                  analysis_data_downsampled$Class_final == "Aves"),
            MAX_PAIRS_PER_GENUS))

genus_dist_after <- analysis_data_downsampled %>%
  count(Class_final, genus_final, name = "n_pairs") %>%
  arrange(Class_final, desc(n_pairs))

cat("\nTop 10 genera by pair count (after downsampling):\n")
print(head(genus_dist_after, 10))
cat(sprintf("Max pairs in any genus: %d (cap = %d)\n",
            max(genus_dist_after$n_pairs), MAX_PAIRS_PER_GENUS))

cat(sprintf("\nPairs before downsampling: %d\n", nrow(analysis_data)))
cat(sprintf("Pairs after downsampling:  %d\n", nrow(analysis_data_downsampled)))

downsampling_report <- genus_dist_before %>%
  left_join(genus_dist_after, by = c("Class_final", "genus_final"),
            suffix = c("_before", "_after")) %>%
  mutate(n_pairs_after  = replace_na(n_pairs_after, 0),
         n_removed      = n_pairs_before - n_pairs_after) %>%
  filter(n_removed > 0) %>%
  arrange(desc(n_removed))

write_csv(downsampling_report, "results/genus_downsampling_report.csv")
cat("\nGenera that were downsampled:\n")
print(downsampling_report)

# use downsampled data from now on
analysis_data <- analysis_data_downsampled

# 12) split by taxonomic group to allow different models

data_by_class <- list(
  mammals    = analysis_data %>% filter(Class_final == "Mammalia"),
  birds      = analysis_data %>% filter(Class_final == "Aves")
)

cat("\n=== DATASET SUMMARIES ===\n")
for (class_name in names(data_by_class)) {
  df <- data_by_class[[class_name]]
  n_pairs <- nrow(df)
  n_species <- length(unique(c(df$sp1, df$sp2)))
  n_genera <- length(unique(df$genus_final))
  cat(sprintf("%s: %d pairs, %d species, %d genera\n", 
              class_name, n_pairs, n_species, n_genera))
}

# check that sample sizes are sufficient for the model complexity
for (class_name in names(data_by_class)) {
  df <- data_by_class[[class_name]]
  n_pairs    <- nrow(df)
  n_species  <- length(unique(c(df$sp1, df$sp2)))
  n_genera   <- length(unique(df$genus_final))
  n_complete <- df %>%
    filter(!is.na(Body.size..kg._avg), !is.na(Generation.time..years._avg),
           !is.na(Clutch.size_avg), !is.na(Population.Density_avg),
           !is.na(Genome.size_avg), !is.na(Genome.GC.content_avg),
           !is.na(N50_bp_avg)) %>%
    nrow()
  cat(sprintf("%s: %d pairs, %d species, %d genera, %d complete cases\n",
              class_name, n_pairs, n_species, n_genera, n_complete))
}

data_by_class <- data_by_class[RUN_CLASS]


# 13) tree processing functions

# normalize tree tip labels to match species names in data.
normalize_tree_labels <- function(tree) {
  tree$tip.label <- normalize_species(tree$tip.label)
  tree
}

# build phylogenetic covariance matrix from pruned tree
build_phylo_cov <- function(tree, chronos_timeout = 300) {
  if (!is.ultrametric(tree)) {
    cat("  Tree not ultrametric, running chronos...\n")
    
    tree <- tryCatch({
      setTimeLimit(elapsed = chronos_timeout, transient = TRUE)
      result <- chronos(tree, quiet = TRUE)
      setTimeLimit(elapsed = Inf, transient = TRUE)
      cat("  chronos completed successfully\n")
      result
    }, error = function(e) {
      setTimeLimit(elapsed = Inf, transient = TRUE)
      if (grepl("time limit", e$message, ignore.case = TRUE)) {
        cat(sprintf("  chronos timed out after %d seconds — trying force.ultrametric\n",
                    chronos_timeout))
      } else {
        cat(sprintf("  chronos failed (%s) — trying force.ultrametric\n", e$message))
      }
      return(NULL)
    })
    
    # fallback if chronos failed or timed out
    if (is.null(tree)) {
      tree <- tryCatch({
        phytools::force.ultrametric(tree, method = "extend")
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
  df_filtered <- df %>%
    filter(sp1 %in% keep_species, sp2 %in% keep_species)
  n_before <- nrow(df)
  n_after  <- nrow(df_filtered)
  cat(sprintf("  Filtered data: %d -> %d pairs (%.1f%% retained)\n",
              n_before, n_after, 100 * n_after / n_before))
  return(df_filtered)
}


# 14) log and z transform predictors

scale_predictors <- function(df) {
  df %>%
    mutate(

      log_k2p = log(k2p),
      log_timetree = log(timetree_div),
      
      # log-transform skewed predictors (add small constant to avoid log(0))
      log_body_size_avg  = log(Body.size..kg._avg + 0.001),
      log_body_size_diff = log(Body.size..kg._diff + 0.001),
      log_gen_time_avg   = log(Generation.time..years._avg + 1),
      log_gen_time_diff  = log(Generation.time..years._diff + 1),
      log_genome_size_avg  = log(Genome.size_avg),
      log_genome_size_diff = log(Genome.size_diff + 0.001),
      log_n50_avg  = log(N50_bp_avg + 1),
      log_n50_diff = log(N50_bp_diff + 1),
      log_area_avg = log(area_km2_avg + 1),
      log_pop_density_avg = log(Population.Density_avg + 1),
      
      # Z-score standardization (mean=0, sd=1)
      z_timetree = scale(log_timetree)[,1],
      z_body_size_avg  = scale(log_body_size_avg)[,1],
      z_body_size_diff = scale(log_body_size_diff)[,1],
      z_gen_time_avg   = scale(log_gen_time_avg)[,1],
      z_gen_time_diff  = scale(log_gen_time_diff)[,1],
      z_clutch_avg     = scale(Clutch.size_avg)[,1],
      z_clutch_diff    = scale(Clutch.size_diff)[,1],
      z_pop_density_avg = scale(log_pop_density_avg)[,1],
      z_genome_size_avg  = scale(log_genome_size_avg)[,1],
      z_genome_size_diff = scale(log_genome_size_diff)[,1],
      z_gc_avg  = scale(Genome.GC.content_avg)[,1],
      z_gc_diff = scale(Genome.GC.content_diff)[,1],
      z_n50_avg  = scale(log_n50_avg)[,1],
      z_n50_diff = scale(log_n50_diff)[,1]
    )
}


# 15) synonym substitution

# synonyms were identified manually by searching the GBIF database 

load_synonyms <- function() {
  bird_syn <- tryCatch(
    read_csv(
      "/scratch/gautschi/allen715/GD_models/bird_synonyms.csv",
      show_col_types = FALSE
    ),
    error = function(e) {
      warning(sprintf("Could not load bird_synonyms.csv: %s", e$message))
      return(NULL)
    }
  )
  
  mammal_syn <- tryCatch(
    read_csv(
      "/scratch/gautschi/allen715/GD_models/mammal_synonyms.csv",
      show_col_types = FALSE
    ),
    error = function(e) {
      warning(sprintf("Could not load mammal_synonyms.csv: %s", e$message))
      return(NULL)
    }
  )
  
  bind_rows(bird_syn, mammal_syn)
}

apply_synonym_substitutions <- function(tree, class_name, dataset_species,
                                        synonyms_df, tree_id, verbose = TRUE) {
  # filter synonyms to this class and to rows that have a synonym in the tree
  class_label <- if (class_name == "birds") "birds" else "mammals"
  
  syn_this_class <- synonyms_df %>%
    filter(
      class == class_label,
      !is.na(`synonym in tree`),
      str_trim(`synonym in tree`) != ""
    ) %>%
    # normalize both columns to lowercase_underscore to match tree tip format
    mutate(
      dataset_name  = normalize_species(`species missing in tree`),
      synonym_norm  = normalize_species(`synonym in tree`)
    )
  
  if (nrow(syn_this_class) == 0) {
    if (verbose) cat("  No synonyms to apply for this class.\n")
    return(list(tree = tree, n_substituted = 0, substitutions = tibble()))
  }
  
  current_tips <- tree$tip.label
  
  # build substitution table
  substitution_log <- syn_this_class %>%
    mutate(
      synonym_in_tree   = synonym_norm %in% current_tips,
      target_in_dataset = dataset_name %in% dataset_species,
      # skip if synonym already appears in the dataset as its own valid species
      synonym_also_in_dataset = synonym_norm %in% dataset_species,
      eligible = synonym_in_tree & !synonym_also_in_dataset
    )
  
  eligible <- substitution_log %>% filter(eligible)
  skipped_conflict  <- substitution_log %>%
    filter(synonym_in_tree & synonym_also_in_dataset)
  skipped_not_found <- substitution_log %>%
    filter(!synonym_in_tree & !is.na(synonym_norm))
  
  # apply substitutions
  n_substituted <- 0
  if (nrow(eligible) > 0) {
    for (i in seq_len(nrow(eligible))) {
      old_name <- eligible$synonym_norm[i]
      new_name <- eligible$dataset_name[i]
      tree$tip.label[tree$tip.label == old_name] <- new_name
      n_substituted <- n_substituted + 1
    }
  }
  
  if (verbose) {
    cat(sprintf("  Synonym substitutions: %d applied", n_substituted))
    if (nrow(skipped_conflict) > 0) {
      cat(sprintf(", %d skipped (synonym also in dataset)", nrow(skipped_conflict)))
    }
    if (nrow(skipped_not_found) > 0) {
      cat(sprintf(", %d synonyms not found in this tree", nrow(skipped_not_found)))
    }
    cat("\n")
    
    if (n_substituted > 0) {
      cat("  Substitutions made:\n")
      for (i in seq_len(nrow(eligible))) {
        cat(sprintf("    %s  →  %s\n",
                    eligible$synonym_norm[i], eligible$dataset_name[i]))
      }
    }
    if (nrow(skipped_conflict) > 0) {
      cat("  Skipped (synonym exists in dataset — would create duplicate):\n")
      for (i in seq_len(nrow(skipped_conflict))) {
        cat(sprintf("    %s (dataset name: %s)\n",
                    skipped_conflict$synonym_norm[i],
                    skipped_conflict$dataset_name[i]))
      }
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

# summarize substitutions made
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
  
  write_csv(sub_freq,
            file.path(out_dir, "synonym_substitution_summary.csv"))
  
  cat("\n=== SYNONYM SUBSTITUTION SUMMARY (ALL TREES) ===\n")
  for (cl in unique(sub_freq$class)) {
    df_cl <- sub_freq %>% filter(class == cl)
    cat(sprintf("\n%s: %d unique substitutions across trees\n",
                toupper(cl), nrow(df_cl)))
    print(df_cl, n = 50)
  }
  
  cat(sprintf("\nFull summary saved to: %s\n",
              file.path(out_dir, "synonym_substitution_summary.csv")))
  
  return(invisible(sub_freq))
}


# 16) phylogenetic signal tests

# calculate Pagel's λ and Blomberg's K to assess phylogenetic structure
# justifies use of phylogenetic mixed models

test_phylogenetic_signal <- function(df, tree, class_name) {
  cat(sprintf("\n=== PHYLOGENETIC SIGNAL: %s ===\n", class_name))
  
  # aggregate k2p values by species (one value per species)
  trait_data <- df %>%
    select(sp1, k2p) %>%
    group_by(sp1) %>%
    summarize(mean_k2p = mean(log(k2p), na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(mean_k2p))
  
  species_in_both <- intersect(tree$tip.label, trait_data$sp1)
  
  if (length(species_in_both) < 10) {
    cat("Insufficient species for phylogenetic signal test\n")
    return(NULL)
  }
  
  tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, species_in_both))
  
  trait_vec <- trait_data %>%
    filter(sp1 %in% species_in_both) %>%
    arrange(match(sp1, tree_pruned$tip.label)) %>%
    pull(mean_k2p)
  names(trait_vec) <- tree_pruned$tip.label
  
  # calculate Pagel's λ (ranges 0-1, tests for phylogenetic signal)
  lambda_result <- tryCatch({
    phylosig(tree_pruned, trait_vec, method = "lambda", test = TRUE)
  }, error = function(e) {
    cat(sprintf("Error calculating lambda: %s\n", e$message))
    return(NULL)
  })
  
  # calculate Blomberg's K (K≈1 means Brownian motion evolution)
  K_result <- tryCatch({
    phylosig(tree_pruned, trait_vec, method = "K", test = TRUE)
  }, error = function(e) {
    cat(sprintf("Error calculating K: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(lambda_result)) {
    cat(sprintf("Pagel's λ = %.3f (p = %.4f)\n", lambda_result$lambda, lambda_result$P))
  }
  
  if (!is.null(K_result)) {
    cat(sprintf("Blomberg's K = %.3f (p = %.4f)\n", K_result$K, K_result$PIC.variance.P))
  }
  
  return(list(lambda = lambda_result, K = K_result))
}

# 17) MODEL FORMULAS BY GROUP

# define predictors for each taxonomic group

get_model_formula <- function(class_name, model_type = "full") {
  
  if (class_name %in% c("mammals", "birds")) {
    # full model with life history and genome traits
    if (model_type == "avg_only") {
      predictors <- "z_timetree + z_body_size_avg + z_gen_time_avg + z_clutch_avg + z_pop_density_avg + z_genome_size_avg + z_gc_avg + z_n50_avg"
    } else if (model_type == "diff_only") {
      predictors <- "z_timetree + z_body_size_diff + z_gen_time_diff + z_clutch_diff + z_genome_size_diff + z_gc_diff + z_n50_diff"
    } else {
      # full model with both average and difference
      predictors <- "z_timetree + z_body_size_avg + z_body_size_diff + z_gen_time_avg + z_clutch_avg + z_pop_density_avg + z_genome_size_avg + z_genome_size_diff + z_gc_avg + z_n50_avg"
    }
  } else {
    stop(sprintf("Unknown taxonomic group: %s", class_name))
  }
  
  # random effects structure:
  # (1 | genus): genus-specific intercepts
  # (1 | mm(sp1, sp2, cov = A)): dyadic phylogenetic effect
  random_effects <- "(1 | genus) + (1 | mm(sp1, sp2, cov = A))"
  
  full_formula <- paste("log_k2p ~", predictors, "+", random_effects)
  
  return(full_formula)
}

# 18) model fitting

# fit Bayesian phylogenetic mixed models using brms/Stan

fit_dyadic_phylo_model <- function(df, A, class_name, model_type = "full") {
  cat(sprintf("\n=== FITTING MODEL: %s (%s) ===\n", class_name, model_type))
  
  # scale predictors
  df <- scale_predictors(df)
  
  # get model formula
  formula_str <- get_model_formula(class_name, model_type)
  cat("\nFormula:\n  ", formula_str, "\n\n")
  
  # extract predictor names to check complete cases
  predictors <- str_extract_all(formula_str, "z_[a-z0-9_]+")[[1]]
  
  # filter to complete cases
  cols_needed <- c("log_k2p", "sp1", "sp2", "genus", predictors)
  df_complete <- df %>%
    select(all_of(cols_needed)) %>%
    na.omit()
  
  n_complete <- nrow(df_complete)
  n_total <- nrow(df)
  pct_complete <- n_complete / n_total
  
  cat(sprintf("Complete cases: %d / %d (%.1f%%)\n", n_complete, n_total, pct_complete * 100))
  
  if (n_complete == 0) {
    stop("No complete cases! Cannot fit model.")
  }
  
  if (n_complete < 20) {
    warning(sprintf("Very few complete cases (%d). Model may not converge.", n_complete))
  }
  
  # fit model using brms (Bayesian Regression Models using Stan)
  cat("\nFitting model (this may take 10-30 minutes)...\n")
  
  fit <- brm(
    formula = as.formula(formula_str),
    data = df_complete,
    data2 = list(A = A),  # phylogenetic covariance matrix
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = b),  # priors for fixed effects
      prior(exponential(1), class = sd),  # priors for random effect SDs
      prior(exponential(1), class = sigma)  # prior for residual SD
    ),
    chains = N_CHAINS,
    cores = N_CHAINS,
    iter = N_ITER,
    warmup = N_WARMUP,
    thin = 2,  # save every 2nd sample to reduce storage
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
    save_pars = save_pars(all = TRUE),  # required for moment matching and LOO
    seed = SEED
  )
  
  return(fit)
}

# 19) model diagnostics

# check convergence and sampling quality

diagnose_model <- function(fit, class_name, model_type = "full") {
  cat(sprintf("\n=== DIAGNOSTICS: %s (%s) ===\n", class_name, model_type))
  
  # Rhat: measures chain convergence (should be < 1.01)
  cat("\nRhat values:\n")
  rhat_vals <- brms::rhat(fit)
  print(summary(rhat_vals))
  
  if (any(rhat_vals > 1.01, na.rm = TRUE)) {
    warning("Some Rhat values > 1.01, model may not have converged!")
  }
  
  # divergent transitions: indicate sampling problems
  cat("\nDivergent transitions:\n")
  div <- brms::nuts_params(fit, pars = "divergent__")
  n_div <- sum(subset(div, Parameter == "divergent__", Value, drop = TRUE))
  cat(sprintf("  %d divergent transitions\n", n_div))
  
  if (n_div > 0) {
    warning("Divergent transitions detected!")
  }
  
  # ESS (Effective Sample Size): measures sampling efficiency
  cat("\nEffective sample size:\n")
  ess <- brms::neff_ratio(fit)
  print(summary(ess))
  
  ess_bulk <- neff_ratio(fit)
  if (any(ess_bulk < 0.1, na.rm = TRUE)) {
    warning("Some parameters have ESS ratio < 0.1")
  }
  
  # return summary diagnostics
  return(list(
    rhat_max = max(rhat_vals, na.rm = TRUE),
    rhat_mean = mean(rhat_vals, na.rm = TRUE),
    ess_min = min(ess, na.rm = TRUE),
    ess_mean = mean(ess, na.rm = TRUE),
    n_divergent = n_div,
    converged = max(rhat_vals, na.rm = TRUE) < 1.01
  ))
}


# 20) sensitivity analysis

# compare models with different random effect structures
# tests whether genus effect is necessary beyond phylogeny

fit_sensitivity_models <- function(df, A, class_name) {
  cat(sprintf("\n=== SENSITIVITY ANALYSIS: %s ===\n", class_name))
  
  # if needed to increase max size
  options(future.globals.maxSize = 3 * 1024^3) 
  
  formula_with_genus <- get_model_formula(class_name, "full")
  formula_no_genus <- str_replace(formula_with_genus, 
                                  "\\+ \\(1 \\| genus\\) \\+ ", "+ ")
  formula_with_family <- str_replace(formula_with_genus,
                                     "\\(1 \\| genus\\)", "(1 | Family)")
  
  df <- scale_predictors(df)
  predictors <- str_extract_all(formula_with_genus, "z_[a-z0-9_]+")[[1]]
  
  cols_needed <- c("log_k2p", "sp1", "sp2", "genus", "Family", predictors)
  df_complete <- df %>%
    select(all_of(cols_needed)) %>%
    na.omit()
  
  cat(sprintf("Fitting sensitivity models with %d complete cases\n", nrow(df_complete)))
  
  # fit Model 1: genus + phylogeny (default)
  cat("\n1. Model with genus effect...\n")
  fit_genus <- brm(
    formula = as.formula(formula_with_genus),
    data = df_complete,
    data2 = list(A = A),
    family = gaussian(),
    prior = c(prior(normal(0, 1), class = b),
              prior(exponential(1), class = sd),
              prior(exponential(1), class = sigma)),
    chains = N_CHAINS, cores = N_CHAINS,
    iter = N_ITER, warmup = N_WARMUP,
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
    seed = SEED
  )
  
  # fit Model 2: phylogeny only
  cat("\n2. Model without genus effect...\n")
  fit_no_genus <- brm(
    formula = as.formula(formula_no_genus),
    data = df_complete,
    data2 = list(A = A),
    family = gaussian(),
    prior = c(prior(normal(0, 1), class = b),
              prior(exponential(1), class = sd),
              prior(exponential(1), class = sigma)),
    chains = N_CHAINS, cores = N_CHAINS,
    iter = N_ITER, warmup = N_WARMUP,
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
    seed = SEED
  )
  
  # fit Model 3: family + phylogeny
  cat("\n3. Model with family effect...\n")
  fit_family <- brm(
    formula = as.formula(formula_with_family),
    data = df_complete,
    data2 = list(A = A),
    family = gaussian(),
    prior = c(prior(normal(0, 1), class = b),
              prior(exponential(1), class = sd),
              prior(exponential(1), class = sigma)),
    chains = N_CHAINS, cores = N_CHAINS,
    iter = N_ITER, warmup = N_WARMUP,
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
    seed = SEED
  )
  
  # compare models using LOO cross-validation
  cat("\n=== MODEL COMPARISON (LOO-CV) ===\n")
  loo_genus    <- loo(fit_genus,    reloo = TRUE, reloo_args = list(chains = 1, cores = 1))
  loo_no_genus <- loo(fit_no_genus, reloo = TRUE, reloo_args = list(chains = 1, cores = 1))
  loo_family   <- loo(fit_family,   reloo = TRUE, reloo_args = list(chains = 1, cores = 1))
  
  comparison <- loo_compare(loo_genus, loo_no_genus, loo_family)
  print(comparison)
  
  cat("\nModel 1 (with genus):\n")
  vc_genus <- VarCorr(fit_genus)
  print(vc_genus)
  
  cat("\nModel 2 (no genus):\n")
  vc_no_genus <- VarCorr(fit_no_genus)
  print(vc_no_genus)
  
  cat("\nModel 3 (with family):\n")
  vc_family <- VarCorr(fit_family)
  print(vc_family)
  
  # save results
  saveRDS(list(
    fit_genus = fit_genus,
    fit_no_genus = fit_no_genus,
    fit_family = fit_family,
    loo_comparison = comparison,
    variance_components = list(
      genus = vc_genus,
      no_genus = vc_no_genus,
      family = vc_family
    )
  ), sprintf("results/sensitivity/%s_random_effects_sensitivity.rds", class_name))
  
  return(list(
    fit_genus = fit_genus,
    fit_no_genus = fit_no_genus,
    fit_family = fit_family,
    comparison = comparison
  ))
}


# 21) load tree selection script
source("tree_selection_utils.R")

# 22) single tree analysis

# process one set of phylogenetic trees (one per taxonomic group)

run_single_tree_analysis <- function(tree_id, tree_config, data_by_class) {
  cat(sprintf("\n\n>>> TREE %d <<<\n", tree_id))
  
  class_names <- names(data_by_class)
  
  # if all classes already checkpointed, skip tree loading entirely
  if (tree_is_complete(tree_id, class_names)) {
    cat(sprintf("  All classes already complete for tree %d — skipping.\n", tree_id))
    result <- list()
    for (cn in class_names) result[[cn]] <- read_checkpoint(tree_id, cn)
    return(result)
  }
  
  # load synonyms once per tree
  synonyms_df <- load_synonyms()
  
  # load trees for this iteration
  cat("\n=== LOADING TREES ===\n")
  loaded_trees <- load_single_trees(
    base_dir   = "/scratch/gautschi/allen715/GD_models/final_dataset/phylosets",
    mammal_n   = tree_config$mammal_tree_number[tree_id],
    bird_file  = tree_config$bird_file[tree_id],
    bird_tree  = tree_config$bird_tree_index[tree_id],
    squam_file = tree_config$squamate_file[tree_id],
    squam_tree = tree_config$squamate_tree_index[tree_id],
    verbose    = TRUE
  )
  
  cat("\n=== LOADED TREES CHECK ===\n")
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
    
    # check if this class is already done for this tree
    existing <- read_checkpoint(tree_id, class_name)
    if (!is.null(existing)) {
      cat(sprintf("  Class %s already complete — loading from checkpoint.\n", class_name))
      results[[class_name]] <- existing
      next
    }
    
    tree_key <- if (class_name == "reptiles") "squamates" else class_name
    tree <- loaded_trees[[tree_key]]
    
    if (is.null(tree)) {
      warning(sprintf("No tree loaded for %s (tree_key: %s)", class_name, tree_key))
      next
    }
    
    cat(sprintf("  Tree loaded: %d tips\n", length(tree$tip.label)))
    
    tree <- normalize_tree_labels(tree)
    
    # apply synonym substitutions before pruning
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
    
    # prune tree to dataset (after synonym substitution)
    pruning_result <- prune_tree_to_dataset(
      tree       = tree,
      df         = data_by_class[[class_name]],
      class_name = class_name,
      tree_id    = tree_id,
      min_tips   = 10,
      verbose    = TRUE
    )
    
    if (is.null(pruning_result$tree)) {
      next
    }
    
    tree_pruned <- pruning_result$tree
    A <- tryCatch(
      build_phylo_cov(tree_pruned),
      error = function(e) {
        cat(sprintf("  build_phylo_cov failed for %s tree %d: %s — skipping\n",
                    class_name, tree_id, e$message))
        return(NULL)
      }
    )

    if (is.null(A)) next

    df <- prepare_brms_data(data_by_class[[class_name]], tree_pruned)
    
    if (nrow(df) < 10) {
      warning(sprintf("Only %d pairs after filtering for %s (tree %d) — skipping.",
                      nrow(df), class_name, tree_id))
      next
    }
    
    cat(sprintf("  Ready to fit model: %d pairs, %d species in tree\n",
                nrow(df), pruning_result$n_tips))
    
    if (tree_id == 1 && RUN_PHYLO_SIGNAL_TESTS) {
      phylo_signal <- test_phylogenetic_signal(df, tree_pruned, class_name)
      saveRDS(phylo_signal, sprintf("results/diagnostics/%s_phylo_signal.rds", class_name))
    }
    
    tryCatch({
      cat("  Fitting average-only model...\n")
      fit_avg  <- fit_dyadic_phylo_model(df, A, class_name, "avg_only")
      diag_avg <- diagnose_model(fit_avg, class_name, "avg_only")
      
      posterior_avg <- as_draws_df(fit_avg)
      saveRDS(posterior_avg,
              sprintf("results/posteriors/%s_tree%04d_avg_only.rds",
                      class_name, tree_id))
      
      cat("  Fitting difference-only model...\n")
      fit_diff  <- fit_dyadic_phylo_model(df, A, class_name, "diff_only")
      diag_diff <- diagnose_model(fit_diff, class_name, "diff_only")
      
      posterior_diff <- as_draws_df(fit_diff)
      saveRDS(posterior_diff,
        sprintf("results/posteriors/%s_tree%04d_diff_only.rds",
                class_name, tree_id))
      
      model_summary <- list(
        avg_model = list(
          formula = formula(fit_avg),
          nobs = nobs(fit_avg),
          diagnostics = diag_avg,
          variance_components = VarCorr(fit_avg, summary = TRUE)
        ),
        diff_model = list(
          formula = formula(fit_diff),
          nobs = nobs(fit_diff),
          diagnostics = diag_diff,
          variance_components = VarCorr(fit_diff, summary = TRUE)
        )
      )
      
      summary_file <- sprintf("results/models/%s_tree%04d_summary.rds",
                              class_name, tree_id)
      saveRDS(model_summary, summary_file)
      cat(sprintf("  ✓ Saved summary: %s\n", summary_file))
      
      class_result <- list(
        avg_model = list(
          posterior   = posterior_avg,
          diagnostics = diag_avg
        ),
        diff_model = list(
          posterior   = posterior_diff,
          diagnostics = diag_diff
        ),
        synonym_log = syn_result
      )
      
      results[[class_name]] <- class_result
      
      # write checkpoint immediately after this class completes
      write_checkpoint(tree_id, class_name, class_result)
      
      cat(sprintf("  ✓ %s model fitted successfully\n", class_name))
      
      # Tree-1-only extras
      if (tree_id == 1 && RUN_SENSITIVITY_ANALYSES) {
        cat("  Running sensitivity analysis...\n")
        fit_sensitivity_models(df, A, class_name)
      }
      
      if (tree_id == 1) {
        saveRDS(
          list(tree = tree_pruned, A = A, df = df),
          sprintf("results/downsampling_sensitivity/%s_tree1_objects.rds", class_name)
        )
      }
      
    }, error = function(e) {
      cat(sprintf("  ✗ Error fitting %s: %s\n", class_name, e$message))
      cat(sprintf("  Full error: %s\n", toString(e)))
      results[[class_name]] <<- list(error = e$message)
    })
  }
  
  return(c(results, list(.synonym_log = synonym_log)))
}


# 23) multi tree pipeline 

# ── TREE SUBSTITUTION POOL ────────────────────────────────────────────────────
# When a tree fails (chronos timeout, pruning failure, etc.), automatically
# draw a replacement from the same posterior distribution so N stays constant.

generate_tree_config_with_substitutes <- function(n_trees, n_reserve = 50, seed = SEED) {
  set.seed(seed)
  # Generate more trees than needed — extras serve as substitutes
  full_config <- generate_tree_sample_config(n_trees = n_trees + n_reserve)
  
  list(
    primary   = full_config[1:n_trees, ],
    reserve   = full_config[(n_trees + 1):(n_trees + n_reserve), ]
  )
}

run_multi_tree_analysis <- function(tree_config_list, data_by_class) {
  
  # accept either old-style single config or new list with reserve pool
  if (is.data.frame(tree_config_list)) {
    tree_config <- tree_config_list
    reserve     <- NULL
  } else {
    tree_config <- tree_config_list$primary
    reserve     <- tree_config_list$reserve
  }
  
  n_trees     <- nrow(tree_config)
  class_names <- names(data_by_class)
  
  cat("\n=================================================\n")
  cat(sprintf("RUNNING ANALYSIS ON %d TREES\n", n_trees))
  if (!is.null(reserve)) {
    cat(sprintf("Reserve pool: %d substitute trees available\n", nrow(reserve)))
  }
  cat(sprintf("Parallelization: %s\n", ifelse(USE_PARALLEL, "ENABLED", "DISABLED")))
  cat("=================================================\n")
  
  n_already_done <- report_checkpoint_status(n_trees, class_names)
  
  if (n_already_done == n_trees) {
    cat("All trees already complete. Consolidating results.\n")
    return(consolidate_checkpoints(n_trees, class_names))
  }
  
  # track which reserve trees have been used
  reserve_used  <- 0
  all_results   <- vector("list", n_trees)
  
  run_one_tree <- function(i, config_row) {
    result <- run_single_tree_analysis(i, config_row, data_by_class)
    
    # check if this tree produced usable results for all classes
    has_error <- sapply(class_names, function(cn) {
      r <- result[[cn]]
      is.null(r) || 
        isTRUE(r$skipped) ||
        (!is.null(r$avg_model) && !is.null(r$avg_model$diagnostics$error))
    })
    
    return(list(result = result, failed_classes = class_names[has_error]))
  }
  
  i <- 1
  reserve_idx <- 1
  
  while (i <= n_trees) {
    
    # skip if already checkpointed
    if (tree_is_complete(i, class_names)) {
      cat(sprintf("  Tree %d already complete — loading from checkpoint.\n", i))
      all_results[[i]] <- consolidate_checkpoints_single(i, class_names)
      i <- i + 1
      next
    }
    
    # build a single-row config for this tree
    config_row <- tree_config[i, , drop = FALSE]
    # inject the tree_id as row index 1 so run_single_tree_analysis addresses it correctly
    rownames(config_row) <- "1"
    
    cat(sprintf("\n>>> TREE %d (primary) <<<\n", i))
    attempt <- run_one_tree(i, rbind(config_row, tree_config))  
    # note: pass full config so tree_id indexing works; run_single_tree_analysis uses tree_id
    
    if (length(attempt$failed_classes) == 0) {
      # success
      all_results[[i]] <- attempt$result
      i <- i + 1
    } else {
      # failure — try a substitute tree
      cat(sprintf("  Tree %d failed for classes: %s\n",
                  i, paste(attempt$failed_classes, collapse = ", ")))
      
      if (!is.null(reserve) && reserve_idx <= nrow(reserve)) {
        cat(sprintf("  Substituting with reserve tree %d (pool index %d)\n",
                    i, reserve_idx))
        
        # replace the failed row in tree_config with the reserve tree
        # keep the same tree_id (i) so checkpoints write to the correct slot
        tree_config[i, ] <- reserve[reserve_idx, ]
        reserve_used      <- reserve_used + 1
        reserve_idx       <- reserve_idx + 1
        
        # log the substitution
        sub_log_path <- "results/checkpoints/substitution_log.csv"
        sub_entry <- data.frame(
          tree_id       = i,
          original_tree = attempt$result$.failed_tree %||% "unknown",
          substitute_from_reserve = reserve_idx - 1,
          timestamp     = as.character(Sys.time())
        )
        if (file.exists(sub_log_path)) {
          write_csv(bind_rows(read_csv(sub_log_path, show_col_types = FALSE),
                              sub_entry), sub_log_path)
        } else {
          write_csv(sub_entry, sub_log_path)
        }
        
        # retry with substitute — do NOT increment i
        cat(sprintf("  Retrying tree slot %d with substitute...\n", i))
        
      } else {
        cat(sprintf("  No reserve trees remaining — tree %d will be missing.\n", i))
        all_results[[i]] <- list(skipped = TRUE, skip_reason = "no_reserve_available")
        i <- i + 1
      }
    }
  }
  
  if (reserve_used > 0) {
    cat(sprintf("\n%d trees were substituted from the reserve pool.\n", reserve_used))
    cat("See results/checkpoints/substitution_log.csv for details.\n")
  }
  
  # synonym and pruning summaries
  all_synonym_logs <- unlist(
    lapply(all_results, function(tree_result) {
      sl <- tree_result[[".synonym_log"]]
      if (is.null(sl)) list() else sl
    }),
    recursive = FALSE
  )
  summarize_synonym_substitutions(all_synonym_logs)
  summarize_tree_pruning(class_names = class_names)
  
  return(all_results)
}

# load a single tree's results from checkpoints
consolidate_checkpoints_single <- function(tree_id, class_names) {
  result <- list()
  for (cn in class_names) {
    cp <- read_checkpoint(tree_id, cn)
    if (!is.null(cp)) result[[cn]] <- cp
  }
  return(result)
}

# null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# 24) resume system

# enables resuming after timeout or failure.
# after each tree completes, a checkpoint file is written to
# results/checkpoints/tree_{id}_{class}.rds
# on resume, completed trees are skipped and results loaded from disk.

CHECKPOINT_DIR <- "results/checkpoints"
dir.create(CHECKPOINT_DIR, showWarnings = FALSE)

checkpoint_path <- function(tree_id, class_name) {
  sprintf("%s/tree%04d_%s.rds", CHECKPOINT_DIR, tree_id, class_name)
}

write_checkpoint <- function(tree_id, class_name, result) {
  path <- checkpoint_path(tree_id, class_name)
  saveRDS(result, path)
  cat(sprintf("  ✓ Checkpoint saved: %s\n", path))
}

read_checkpoint <- function(tree_id, class_name) {
  path <- checkpoint_path(tree_id, class_name)
  if (file.exists(path)) {
    cat(sprintf("  ↩ Resuming from checkpoint: %s\n", path))
    return(readRDS(path))
  }
  return(NULL)
}

tree_is_complete <- function(tree_id, class_names) {
  all(sapply(class_names, function(cn) file.exists(checkpoint_path(tree_id, cn))))
}

# scan checkpoints and report resume status at startup
report_checkpoint_status <- function(n_trees, class_names) {
  completed <- sum(sapply(1:n_trees, tree_is_complete, class_names = class_names))
  if (completed == 0) {
    cat("  No checkpoints found — starting fresh.\n")
  } else {
    cat(sprintf("  Found checkpoints for %d / %d trees — will resume from tree %d.\n",
                completed, n_trees, completed + 1))
    # list any partially completed trees (some classes done, not all)
    partial <- sapply(1:n_trees, function(i) {
      done <- sapply(class_names, function(cn) file.exists(checkpoint_path(i, cn)))
      any(done) && !all(done)
    })
    if (any(partial)) {
      cat(sprintf("  Partially completed trees (will re-run missing classes): %s\n",
                  paste(which(partial), collapse = ", ")))
    }
  }
  return(completed)
}

# consolidate all checkpoint results into a single list
# called at the end of run_multi_tree_analysis() to assemble final results
consolidate_checkpoints <- function(n_trees, class_names) {
  cat("\n=== CONSOLIDATING RESULTS FROM CHECKPOINTS ===\n")
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

# 24) run models

# set seed before generating config so tree selection is identical if the job is restarted (same trees, same order)
set.seed(SEED)
tree_config_path      <- sprintf("results/tree_config_n%d.csv", N_TREES)
tree_reserve_path     <- sprintf("results/tree_config_n%d_reserve.csv", N_TREES)

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

save.image(file = sprintf("results/workspace_n%d_%s_trees.RData", N_TREES, RUN_CLASS))

