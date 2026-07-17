library(ape)
library(phytools)   # force.ultrametric
library(dplyr)

RES <- "/scratch/gautschi/allen715/GD_models/rate_model/results"   # your results dir
dir.create(file.path(RES, "downsampling_sensitivity"),
           showWarnings = FALSE, recursive = TRUE)

# rebuild the phylogenetic correlation matrix from a (possibly non-ultrametric) tree
build_A <- function(tree) {
  if (!is.ultrametric(tree)) tree <- force.ultrametric(tree, method = "extend")
  vcv.phylo(tree, corr = TRUE)
}

for (cn in c("birds", "mammals")) {
  tree_path <- file.path(RES, "tree_objects", sprintf("%s_tree0001.rds", cn))
  sdf_path  <- file.path(RES, "models",       sprintf("%s_scaled_df.rds", cn))
  if (!file.exists(tree_path) || !file.exists(sdf_path)) {
    cat(sprintf("!! %s: missing %s or %s — skipping\n", cn, tree_path, sdf_path)); next
  }
  tree <- readRDS(tree_path)
  sdf  <- readRDS(sdf_path)
  # keep only pairs whose two species are both in this tree (= prepare_brms_data)
  df   <- sdf %>% filter(sp1 %in% tree$tip.label, sp2 %in% tree$tip.label)
  out  <- file.path(RES, "downsampling_sensitivity", sprintf("%s_tree1_objects.rds", cn))
  saveRDS(list(tree = tree, A = build_A(tree), df = df), out)
  cat(sprintf("wrote %s  (%d tips, %d pairs)\n", out, length(tree$tip.label), nrow(df)))
}
cat("\nDone. Now rerun variance_partitioning.R as-is.\n")
