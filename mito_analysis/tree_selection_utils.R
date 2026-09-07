### tree selection utilities used by pmm scripts to load phylosets

library(ape)

load_single_trees <- function(
    base_dir = "/scratch/gautschi/allen715/GD_models/final_dataset/phylosets",
    mammal_n = NULL,
    bird_file = NULL,
    bird_tree = NULL,
    squam_file = NULL,
    squam_tree = NULL,
    amph_file = NULL,
    amph_tree = NULL,
    verbose = TRUE
) {
  
  trees <- list()
  
  ### mammals
  
  if (!is.null(mammal_n)) {
    if (verbose) cat("\n>>> Loading MAMMAL tree...\n")
    
    mammal_dir <- file.path(base_dir, "DNAonly_4098sp_topoFree_FBDasZhouEtAl")
    
    if (!dir.exists(mammal_dir)) {
      warning(sprintf("Mammal directory not found: %s", mammal_dir))
      trees$mammals <- NULL
    } else {
      mammal_files <- list.files(mammal_dir, pattern = ".*\\.tre$", full.names = TRUE)
      
      if (length(mammal_files) == 0) {
        warning(sprintf("No mammal tree files found in %s", mammal_dir))
        trees$mammals <- NULL
      } else {
        mammal_files <- sort(mammal_files)
        
        if (mammal_n > length(mammal_files)) {
          warning(sprintf("Mammal tree %d exceeds available files (%d)", 
                          mammal_n, length(mammal_files)))
          trees$mammals <- NULL
        } else {
          file_path <- mammal_files[mammal_n]
          if (verbose) cat(sprintf("  File: %s\n", basename(file_path)))
          
          tree <- tryCatch({
            read.tree(file_path)
          }, error = function(e) {
            if (verbose) cat("  read.tree failed, trying read.nexus...\n")
            read.nexus(file_path)
          })
          
          # handle multiPhylo 
          if (inherits(tree, "multiPhylo")) {
            tree <- tree[[1]]
            if (verbose) cat("  Extracted first tree from multiPhylo\n")
          }
          
          if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree$tip.label)))
          trees$mammals <- tree
        }
      }
    }
  }
  
  ### birds
  
  if (!is.null(bird_file) && !is.null(bird_tree)) {
    if (verbose) cat("\n>>> Loading BIRD tree...\n")
    
    # determine path
    if (bird_file == 1) {
      bird_path <- file.path(base_dir, "mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1.tre")
      if (!file.exists(bird_path)) {
        bird_path <- file.path(base_dir, "mnt-1/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1.tre")
      }
    } else {
      bird_path <- file.path(
        base_dir,
        sprintf("mnt-%d/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_%d.tre",
                bird_file, bird_file)
      )
    }
    
    if (!file.exists(bird_path)) {
      warning(sprintf("Bird tree file not found: %s", bird_path))
      trees$birds <- NULL
    } else {
      if (verbose) cat(sprintf("  File: %s\n", basename(bird_path)))
      
      tree_data <- tryCatch({
        read.tree(bird_path)
      }, error = function(e) {
        if (verbose) cat("  read.tree failed, trying read.nexus...\n")
        read.nexus(bird_path)
      })
      
      # extract tree from multiPhylo
      if (inherits(tree_data, "multiPhylo")) {
        n_trees <- length(tree_data)
        if (verbose) cat(sprintf("  File contains %d trees\n", n_trees))
        
        if (bird_tree > n_trees) {
          warning(sprintf("Bird tree index %d exceeds available trees (%d)", 
                          bird_tree, n_trees))
          trees$birds <- NULL
        } else {
          tree <- tree_data[[bird_tree]]
          if (verbose) cat(sprintf("  Using tree #%d\n", bird_tree))
          if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree$tip.label)))
          trees$birds <- tree
        }
      } else {
        if (verbose) cat("  File contains 1 tree\n")
        if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree_data$tip.label)))
        trees$birds <- tree_data
      }
    }
  }
  
  ### reptiles (not used in manuscript)
  
  if (!is.null(squam_file) && !is.null(squam_tree)) {
    if (verbose) cat("\n>>> Loading SQUAMATE tree...\n")
    
    squamate_path <- file.path(
      base_dir,
      sprintf("squam_shl_new_Posterior_9755/squam_shl_new_Posterior_9755.%d000.trees",
              squam_file)
    )
    
    if (!file.exists(squamate_path)) {
      warning(sprintf("Squamate tree file not found: %s", squamate_path))
      trees$squamates <- NULL
    } else {
      if (verbose) cat(sprintf("  File: %s\n", basename(squamate_path)))
      
      tree_data <- tryCatch({
        read.nexus(squamate_path)
      }, error = function(e) {
        if (verbose) cat("  read.nexus failed, trying read.tree...\n")
        read.tree(squamate_path)
      })
      
      # extract tree
      if (inherits(tree_data, "multiPhylo")) {
        n_trees <- length(tree_data)
        if (verbose) cat(sprintf("  File contains %d trees\n", n_trees))
        
        if (squam_tree > n_trees) {
          warning(sprintf("Squamate tree index %d exceeds available trees (%d)", 
                          squam_tree, n_trees))
          trees$squamates <- NULL
        } else {
          tree <- tree_data[[squam_tree]]
          if (verbose) cat(sprintf("  Using tree #%d\n", squam_tree))
          if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree$tip.label)))
          trees$squamates <- tree
        }
      } else {
        if (verbose) cat("  File contains 1 tree\n")
        if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree_data$tip.label)))
        trees$squamates <- tree_data
      }
    }
  }
  
  ### amphibians (not used in manuscript)
  
  if (!is.null(amph_file) && !is.null(amph_tree)) {
    if (verbose) cat("\n>>> Loading AMPHIBIAN tree...\n")
    
    amphibian_path <- file.path(
      base_dir,
      sprintf("amph_shl_new_Posterior_7238/amph_shl_new_Posterior_7238.%d000.trees",
              amph_file)
    )
    
    if (!file.exists(amphibian_path)) {
      warning(sprintf("Amphibian tree file not found: %s", amphibian_path))
      trees$amphibians <- NULL
    } else {
      if (verbose) cat(sprintf("  File: %s\n", basename(amphibian_path)))
      
      tree_data <- tryCatch({
        read.nexus(amphibian_path)
      }, error = function(e) {
        if (verbose) cat("  read.nexus failed, trying read.tree...\n")
        read.tree(amphibian_path)
      })
      
      # extract tree
      if (inherits(tree_data, "multiPhylo")) {
        n_trees <- length(tree_data)
        if (verbose) cat(sprintf("  File contains %d trees\n", n_trees))
        
        if (amph_tree > n_trees) {
          warning(sprintf("Amphibian tree index %d exceeds available trees (%d)", 
                          amph_tree, n_trees))
          trees$amphibians <- NULL
        } else {
          tree <- tree_data[[amph_tree]]
          if (verbose) cat(sprintf("  Using tree #%d\n", amph_tree))
          if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree$tip.label)))
          trees$amphibians <- tree
        }
      } else {
        if (verbose) cat("  File contains 1 tree\n")
        if (verbose) cat(sprintf("  Tree has %d tips\n", length(tree_data$tip.label)))
        trees$amphibians <- tree_data
      }
    }
  }
  
  if (verbose) cat("\n")
  return(trees)
}

### generate tree config

generate_tree_sample_config <- function(n_trees, seed = 12345, base_dir = "/scratch/gautschi/allen715/GD_models/final_dataset/phylosets") {
  set.seed(seed)
  
  cat(sprintf("\nGenerating configuration for %d trees per taxonomic group...\n", n_trees))
  
  # get available mammal tree files
  mammal_dir <- file.path(base_dir, "DNAonly_4098sp_topoFree_FBDasZhouEtAl")
  
  if (!dir.exists(mammal_dir)) {
    warning(sprintf("Mammal directory not found: %s", mammal_dir))
    n_mammal_files <- 10000  # assume max
  } else {
    mammal_files <- list.files(mammal_dir, pattern = "\\.tre$", full.names = FALSE)
    n_mammal_files <- length(mammal_files)
    cat(sprintf("Found %d mammal tree files\n", n_mammal_files))
  }
  
  config <- data.frame(
    tree_id = 1:n_trees,
    
    mammal_tree_number = if (n_trees <= n_mammal_files) {
      sample(1:n_mammal_files, n_trees, replace = FALSE)
    } else {
      sample(1:n_mammal_files, n_trees, replace = TRUE)
    },
    
    bird_file = sample(1:10, n_trees, replace = TRUE),
    bird_tree_index = sample(1:1000, n_trees, replace = TRUE),
    
    squamate_file = sample(1:10, n_trees, replace = TRUE),
    squamate_tree_index = sample(1:1000, n_trees, replace = TRUE),
    
    amphibian_file = sample(1:10, n_trees, replace = TRUE),
    amphibian_tree_index = sample(1:1000, n_trees, replace = TRUE)
  )
  
  write.csv(config, sprintf("/scratch/gautschi/allen715/GD_models/results/tree_config_n%d.csv", n_trees), row.names = FALSE)
  cat(sprintf("Configuration saved to: /scratch/gautschi/allen715/GD_models/results/tree_config_n%d.csv\n", n_trees))
  
  return(config)
}

### load by row

load_trees_by_config_row <- function(config_row, base_dir = "/scratch/gautschi/allen715/GD_models/final_dataset/phylosets", verbose = TRUE) {
  
  load_single_trees(
    base_dir = base_dir,
    mammal_n = config_row$mammal_tree_number,
    bird_file = config_row$bird_file,
    bird_tree = config_row$bird_tree_index,
    squam_file = config_row$squamate_file,
    squam_tree = config_row$squamate_tree_index,
    amph_file = if ("amphibian_file" %in% names(config_row)) config_row$amphibian_file else NULL,
    amph_tree = if ("amphibian_tree_index" %in% names(config_row)) config_row$amphibian_tree_index else NULL,
    verbose = verbose
  )
}

### testing function

test_tree_loading <- function(base_dir = "/scratch/gautschi/allen715/GD_models/final_dataset/phylosets") {
  
  # test loading one tree from each group
  test_trees <- load_single_trees(
    base_dir = base_dir,
    mammal_n = 1,
    bird_file = 1,
    bird_tree = 1,
    squam_file = 1,
    squam_tree = 1,
    verbose = TRUE
  )
  
  for (group in names(test_trees)) {
    if (is.null(test_trees[[group]])) {
      cat(sprintf("%s: FAILED (NULL)\n", group))
    } else {
      cat(sprintf("%s: SUCCESS (%d tips)\n", group, length(test_trees[[group]]$tip.label)))
    }
  }
  
  return(test_trees)
}

### prune tree to dataset species

prune_tree_to_dataset <- function(tree, df, class_name, tree_id,
                                  min_tips = 10, verbose = TRUE) {
  
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  
  # set of species in filtered dataset
  species_in_data <- unique(c(df$sp1, df$sp2))
  n_dataset <- length(species_in_data)
  
  # species in dataset not in this tree
  not_in_tree <- setdiff(species_in_data, tree$tip.label)
  
  # species in tree that aren't needed — drop 
  extra_in_tree <- setdiff(tree$tip.label, species_in_data)
  tree_pruned <- drop.tip(tree, extra_in_tree)
  
  # species in pruned tree
  in_tree <- tree_pruned$tip.label
  n_tips  <- length(in_tree)
  
  if (verbose) {
    cat(sprintf("  [%s tree %d] Dataset: %d species | In tree: %d | Missing from tree: %d\n",
                class_name, tree_id, n_dataset, n_tips, length(not_in_tree)))
    if (length(not_in_tree) > 0 && verbose) {
      missing_preview <- paste(head(not_in_tree, 5), collapse = ", ")
      if (length(not_in_tree) > 5) missing_preview <- sprintf("%s ... (+%d more)", missing_preview, length(not_in_tree) - 5)
      cat(sprintf("    Missing: %s\n", missing_preview))
    }
  }
  
  ### per-tree species coverage report
  
  out_dir <- "/scratch/gautschi/allen715/GD_models/results/species_lists"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_csv(
    data.frame(
      tree_id    = tree_id,
      class      = class_name,
      species    = sort(in_tree),
      in_dataset = in_tree %in% species_in_data  # should always be TRUE after pruning
    ),
    sprintf("%s/%s_tree%04d_in_tree.csv", out_dir, class_name, tree_id)
  )
  
  if (length(not_in_tree) > 0) {
    write_csv(
      data.frame(
        tree_id = tree_id,
        class   = class_name,
        species = sort(not_in_tree)
      ),
      sprintf("%s/%s_tree%04d_missing_from_tree.csv", out_dir, class_name, tree_id)
    )
  }
  
  summary_file <- sprintf("%s/%s_pruning_summary_all_trees.csv", out_dir, class_name)
  
  summary_row <- data.frame(
    tree_id          = tree_id,
    class            = class_name,
    n_dataset_species = n_dataset,
    n_in_tree        = n_tips,
    n_missing        = length(not_in_tree),
    pct_coverage     = round(100 * n_tips / n_dataset, 1),
    missing_species  = if (length(not_in_tree) == 0) NA_character_
    else paste(sort(not_in_tree), collapse = "|")
  )
  
  if (file.exists(summary_file)) {
    write_csv(summary_row, summary_file, append = TRUE)
  } else {
    write_csv(summary_row, summary_file)
  }
  
  # return NULL tree if below minimum tips threshold
  if (n_tips < min_tips) {
    warning(sprintf("[%s tree %d] Only %d tips after pruning — skipping model fit.",
                    class_name, tree_id, n_tips))
    return(list(tree = NULL, n_tips = n_tips, in_tree = in_tree,
                not_in_tree = not_in_tree, n_dataset = n_dataset))
  }
  
  return(list(
    tree        = tree_pruned,
    n_tips      = n_tips,
    in_tree     = in_tree,
    not_in_tree = not_in_tree,
    n_dataset   = n_dataset
  ))
}

summarize_tree_pruning <- function(class_names = c("mammals", "birds"),
                                   out_dir = "/scratch/gautschi/allen715/GD_models/results/species_lists") {
    
  for (class_name in class_names) {
    summary_file <- sprintf("%s/%s_pruning_summary_all_trees.csv", out_dir, class_name)
    
    if (!file.exists(summary_file)) {
      cat(sprintf("  No pruning summary found for %s\n", class_name))
      next
    }
    
    summary_df <- read_csv(summary_file, show_col_types = FALSE)
    
    cat(sprintf("\n--- %s (%d trees) ---\n", toupper(class_name), nrow(summary_df)))
    cat(sprintf("  Mean coverage: %.1f%% of dataset species per tree\n",
                mean(summary_df$pct_coverage, na.rm = TRUE)))
    cat(sprintf("  Coverage range: %.1f%% – %.1f%%\n",
                min(summary_df$pct_coverage), max(summary_df$pct_coverage)))
    cat(sprintf("  Mean missing per tree: %.1f species\n",
                mean(summary_df$n_missing, na.rm = TRUE)))
    
    all_missing <- summary_df %>%
      filter(!is.na(missing_species)) %>%
      pull(missing_species) %>%
      strsplit("\\|") %>%
      unlist()
    
    if (length(all_missing) > 0) {
      missing_freq <- sort(table(all_missing), decreasing = TRUE)
      n_trees_total <- nrow(summary_df)
      
      cat(sprintf("\n  Species missing from >50%% of trees:\n"))
      frequent_missing <- missing_freq[missing_freq > 0.5 * n_trees_total]
      if (length(frequent_missing) > 0) {
        for (sp in names(frequent_missing)) {
          cat(sprintf("    %s: missing from %d/%d trees (%.0f%%)\n",
                      sp, frequent_missing[[sp]], n_trees_total,
                      100 * frequent_missing[[sp]] / n_trees_total))
        }
      } else {
        cat("    None — all species present in >50% of trees\n")
      }
      
      missing_freq_df <- data.frame(
        species       = names(missing_freq),
        n_trees_missing = as.integer(missing_freq),
        n_trees_total   = n_trees_total,
        pct_trees_missing = round(100 * as.integer(missing_freq) / n_trees_total, 1)
      ) %>% arrange(desc(n_trees_missing))
      
      write_csv(missing_freq_df,
                sprintf("%s/%s_species_missing_frequency.csv", out_dir, class_name))
      cat(sprintf("\n  Full missing frequency saved to: %s_species_missing_frequency.csv\n",
                  class_name))
    } else {
      cat("  No species were missing from any tree.\n")
    }
  }
}