# color palette and settings for manuscript figures

suppressPackageStartupMessages({
  library(tidyverse) 
  library(readxl)
  library(patchwork)  
  library(ggrepel)  
})

### paths 
BASE      <- "."  # project root
RAW_DIV   <- file.path(BASE, "Complete_divergence_05-14-26.csv")
TRAITS    <- file.path(BASE, "Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
YARD      <- file.path(BASE, "yardstick_results")  # yardstick results
PMM       <- file.path(BASE, "PMM")  # PMM results
FIG_DIR   <- file.path(BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

NUC_PAL <- c(
  Mammalia = "#0072B2", 
  Aves     = "#D55E00",  
  Reptilia = "#009E73", 
  Amphibia = "#CC79A7"  
)

MITO_PAL <- c(
  Mammalia = "#7BBDE3",  
  Aves     = "#F2A66F", 
  Reptilia = "#75C9AD",  
  Amphibia = "#E2A9C8"  
)

CLASS_PAL <- NUC_PAL

CC_PAL <- c(
  Mammalia.Nuclear = NUC_PAL[["Mammalia"]], Mammalia.Mitochondrial = MITO_PAL[["Mammalia"]],
  Aves.Nuclear     = NUC_PAL[["Aves"]],     Aves.Mitochondrial     = MITO_PAL[["Aves"]],
  Reptilia.Nuclear = NUC_PAL[["Reptilia"]], Reptilia.Mitochondrial = MITO_PAL[["Reptilia"]],
  Amphibia.Nuclear = NUC_PAL[["Amphibia"]], Amphibia.Mitochondrial = MITO_PAL[["Amphibia"]]
)

LAB_NUC  <- setNames(unname(NUC_PAL[c("Aves","Mammalia","Reptilia","Amphibia")]),
                     c("Birds","Mammals","Reptiles","Amphibians"))
LAB_MITO <- setNames(unname(MITO_PAL[c("Aves","Mammalia","Reptilia","Amphibia")]),
                     c("Birds","Mammals","Reptiles","Amphibians"))

# accent colors
ACCENT <- c(yellow = "#F0E442", blue = "#0072B2",
            pink   = "#CC79A7", green = "#E69F00", grey = "#8A8A8A")

# outlier colors
OUT_COL <- c(`Higher than expected` = "#E69F00", `Lower than expected` = "#6A1B9A")
CASE_FILL <- "#D70040"  # case studies

# variance-partition colors
COMP_PAL <- c(
  `Pure traits`      = "#0072B2",
  `Shared`           = "#CC79A7",
  `Pure phylogeny`   = "#E69F00",
  `Residual`         = "#BDBDBD"
)

CLASS_LABELLER <- c(Aves = "Birds", Mammalia = "Mammals",
                    Reptilia = "Reptiles", Amphibia = "Amphibians",
                    birds = "Birds", mammals = "Mammals")

theme_pub <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(face = "bold", size = rel(1.15), hjust = 0),
      plot.subtitle   = element_text(size = rel(0.9), colour = "grey35",
                                      margin = margin(b = 6)),
      plot.title.position = "plot",
      axis.title      = element_text(size = rel(0.95), colour = "grey15"),
      axis.text       = element_text(size = rel(0.82), colour = "grey25"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.35),
      panel.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      legend.title     = element_text(size = rel(0.85)),
      legend.key.height = unit(10, "pt"),
      strip.text       = element_text(face = "bold", size = rel(0.9),
                                      colour = "grey15"),
      plot.margin      = margin(10, 12, 10, 10)
    )
}
theme_set(theme_pub())

save_fig <- function(plot, stem, width, height, dpi = 350) {
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = pdf_dev)
  tryCatch(
    ggsave(file.path(FIG_DIR, paste0(stem, ".png")), plot,
           width = width, height = height, dpi = dpi,
           bg = "white",
           device = if (capabilities("cairo")) "png" else "png"),
    error = function(e)
      message(sprintf("  PNG skipped for %s (%s)", stem, conditionMessage(e)))
  )
  message(sprintf("  saved %s (%.1f x %.1f in)", stem, width, height))
  invisible(plot)
}

MIN_QUERY_ALIGN <- 80
MIN_REF_ALIGN   <- 80
FOUR_CLASSES    <- c("Aves", "Mammalia", "Reptilia", "Amphibia")

normalize_species <- function(x) {
  x %>% str_trim() %>% str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>% str_replace_all("[^A-Za-z ]", "") %>%
    str_to_lower() %>% str_replace_all(" ", "_")
}

build_pairs <- function() {
  traits <- read_excel(TRAITS)
  colnames(traits) <- make.names(colnames(traits))
  div <- read_csv(RAW_DIV, show_col_types = FALSE)
  colnames(div) <- make.names(colnames(div))

  excl <- traits %>% filter(Exclude. %in% c("yes", "Yes")) %>% pull(accession)

  taxo <- traits %>%
    dplyr::select(accession, Class, Order, Family, genus) %>%
    distinct(accession, .keep_all = TRUE)

  pairs <- div %>%
    filter(!genome1 %in% excl, !genome2 %in% excl,
           !is.na(query_alignment_percent), !is.na(ref_alignment_percent),
           query_alignment_percent >= MIN_QUERY_ALIGN,
           ref_alignment_percent   >= MIN_REF_ALIGN) %>%
    left_join(taxo %>% rename_with(~ paste0(.x, "_1")),
              by = c("genome1" = "accession_1")) %>%
    left_join(taxo %>% rename_with(~ paste0(.x, "_2")),
              by = c("genome2" = "accession_2")) %>%
    mutate(
      Class_final  = coalesce(Class_1,  Class_2),
      Order_final  = coalesce(Order_1,  Order_2),
      Family_final = coalesce(Family_1, Family_2),
      genus_final  = str_to_lower(coalesce(genus_1, genus_2)),
      sp1 = normalize_species(species1), sp2 = normalize_species(species2)
    ) %>%
    filter(Class_final %in% FOUR_CLASSES, !is.na(k2p), !is.na(timetree_div),
           k2p > 0, timetree_div > 0) %>%
    mutate(log_k2p = log(k2p), log_div = log(timetree_div),
           Class = factor(Class_final, levels = FOUR_CLASSES))

  ct <- read_csv(file.path(YARD, "class_yardstick_table.csv"),
                 show_col_types = FALSE) %>%
    dplyr::select(class, intercept, slope, residual_sd)

  pairs <- pairs %>%
    left_join(ct, by = c("Class_final" = "class")) %>%
    mutate(
      expected_log = intercept + slope * log_div,
      residual     = log_k2p - expected_log,
      std_residual = residual / residual_sd
    )
  pairs
}

add_class_label <- function(df, col = "Class_final") {
  df[["ClassLab"]] <- factor(unname(CLASS_LABELLER[as.character(df[[col]])]),
                             levels = c("Birds", "Mammals", "Reptiles", "Amphibians"))
  df
}
LAB_PAL <- LAB_NUC  

MITO   <- file.path(BASE, "mito_pipeline/results")  # mitochondrial result CSVs
MEMBER <- file.path(BASE, "dataset_membership.csv")  # nuclear/mito pair membership

build_mito_pairs <- function() {
  ct <- read_csv(file.path(MITO, "class_yardstick_table.csv"), show_col_types = FALSE) %>%
    dplyr::select(class, intercept, slope, residual_sd)
  scored <- file.path(MITO, "yardstick_pairs_all.csv")
  if (file.exists(scored)) {
    df <- read_csv(scored, show_col_types = FALSE) %>%
      mutate(Class_final = as.character(Class_final),
             genus_final = str_to_lower(genus_final),
             k2p         = as.numeric(k2p),
             log_k2p     = if ("log_k2p" %in% names(.)) log_k2p else log(k2p),
             log_div     = if ("log_div_time" %in% names(.)) log_div_time else log(timetree_div))
    return(df %>%
      left_join(ct, by = c("Class_final" = "class")) %>%
      mutate(expected_log = intercept + slope * log_div,
             residual     = log_k2p - expected_log,
             std_residual = residual / residual_sd,
             Class        = factor(Class_final, levels = FOUR_CLASSES)))
  }
  warning("yardstick_pairs_all.csv not found in ", MITO,
          " -- using approximate reconstruction (export it from mito_yardstick.R for an exact match).")
  cmp <- read_csv(file.path(MITO, "mito_k2p_compare.csv"), show_col_types = FALSE) %>%
    filter(keep == 1, !is.na(k2p_coding), k2p_coding > 0)
  mem <- read_csv(MEMBER, show_col_types = FALSE) %>%
    dplyr::select(pair_id, timetree_div, in_mito_yardstick)
  gmap <- read_excel(TRAITS)
  colnames(gmap) <- make.names(colnames(gmap))
  gmap <- gmap %>% transmute(genus = str_to_lower(genus), Order, Family) %>%
    filter(!is.na(genus)) %>% distinct(genus, .keep_all = TRUE)
  cmp %>%
    left_join(mem, by = "pair_id") %>%
    filter(in_mito_yardstick == 1, !is.na(timetree_div), timetree_div > 0) %>%
    mutate(Class_final = taxonomic_group,
           genus_final = str_to_lower(word(species1, 1)),
           sp1 = normalize_species(species1), sp2 = normalize_species(species2),
           k2p = k2p_coding, log_k2p = log(k2p_coding), log_div = log(timetree_div)) %>%
    left_join(gmap, by = c("genus_final" = "genus")) %>%
    rename(Order_final = Order, Family_final = Family) %>%
    left_join(ct, by = c("Class_final" = "class")) %>%
    mutate(expected_log = intercept + slope * log_div,
           residual     = log_k2p - expected_log,
           std_residual = residual / residual_sd,
           Class        = factor(Class_final, levels = FOUR_CLASSES))
}

make_grid <- function(ct, classes = c("Aves", "Mammalia", "Reptilia"),
                      pal = LAB_PAL) {
  ct %>%
    filter(class %in% classes) %>%
    rowwise() %>%
    do({
      r <- .
      lt <- seq(log(r$div_time_min), log(r$div_time_max), length.out = 200)
      tibble(class = r$class, log_div = lt,
             fit = r$intercept + r$slope * lt,
             lwr = fit - 1.96 * r$residual_sd,
             upr = fit + 1.96 * r$residual_sd)
    }) %>% ungroup() %>%
    mutate(ClassLab = factor(unname(CLASS_LABELLER[class]), levels = names(pal)))
}

dodge_forest <- function(df, yvar, groupvar, width = 0.6) {
  yf <- factor(df[[yvar]]); gf <- factor(df[[groupvar]])
  ng <- nlevels(gf)
  off <- ((seq_len(ng) - 1) - (ng - 1) / 2) * (width / ng)
  names(off) <- levels(gf)
  df$.y <- as.integer(yf) + off[as.character(gf)]
  list(df = df, breaks = seq_len(nlevels(yf)), labels = levels(yf))
}
