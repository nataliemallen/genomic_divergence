# supplementary mitochondrial figures
source("figures_setup.R")

FIG_DIR <- file.path(BASE, "figures_si_mito"); dir.create(FIG_DIR, showWarnings = FALSE)
rd <- function(p) readr::read_csv(p, show_col_types = FALSE)
COMP_PAL2 <- c(`Pure traits` = "#FF7D00", `Shared` = "#15616D",
               `Pure phylogeny` = "#E7DFC6", `Residual` = "#392C3A")
MEMBER_PAL <- c(`Genome only` = "#3B3B3B", `Both` = "#CC79A7", `Mito only` = "#C0392B")
OUT_COL   <- c(`Higher than expected` = "#E69F00", `Lower than expected` = "#6A1B9A")  # dark purple

SHAPE_CRED <- scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                                 labels = c(`TRUE` = "credible", `FALSE` = "not credible"),
                                 name = NULL)

mito <- build_mito_pairs() %>% add_class_label()
XB <- log(c(0.1, 0.3, 1, 3, 10, 30)); XL <- c("0.1","0.3","1","3","10","30")
YBM <- log(c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)); YLM <- c("0.001","0.003","0.01","0.03","0.1","0.3")

ct_m   <- rd(file.path(MITO, "class_yardstick_table.csv"))
grid_m <- make_grid(ct_m, pal = LAB_MITO)
flagged <- mito %>%
  mutate(flag = case_when(std_residual >  1.96 ~ "Higher than expected",
                          std_residual < -1.96 ~ "Lower than expected",
                          TRUE ~ "On the yardstick"))
pick <- function(g1, g2) mito %>% filter((grepl(g1, sp1) & grepl(g2, sp2)) |
                                         (grepl(g2, sp1) & grepl(g1, sp2)))

cases <- bind_rows(
  pick("alouatta_juara", "alouatta_macconnelli")         %>% mutate(label = "Alouatta juara / A. macconnelli (concordant excess)"),
  pick("darevskia_armeniaca", "darevskia_mixta")         %>% mutate(label = "Darevskia armeniaca / D. mixta (mito capture)"),
  pick("papio_anubis", "papio_cynocephalus")             %>% mutate(label = "Papio anubis / P. cynocephalus (introgression)"),
  pick("odocoileus_hemionus", "odocoileus_virginianus")  %>% mutate(label = "Odocoileus hemionus / O. virginianus (introgression)")
) %>% distinct()
figS10_base <- ggplot() +
  geom_line(data = grid_m, aes(log_div, fit, group = ClassLab), colour = "grey70", linewidth = 0.7) +
  geom_point(data = filter(flagged, flag == "On the yardstick"),
             aes(log_div, log_k2p), colour = "grey80", size = 0.8, alpha = 0.5, stroke = 0) +
  geom_point(data = filter(flagged, flag != "On the yardstick"),
             aes(log_div, log_k2p, colour = flag), size = 2.0, alpha = 0.9, stroke = 0) +
  scale_colour_manual(values = OUT_COL, name = NULL) +
  geom_point(data = cases, aes(log_div, log_k2p), shape = 21, fill = "#D70040",
             colour = "grey15", size = 3.4, stroke = 0.7) +
  scale_x_continuous(breaks = XB, labels = XL) +
  scale_y_continuous(breaks = YBM, labels = YLM) +
  labs(title = "The mitochondrial yardstick flags introgression and deep splits",
       subtitle = "Pairs outside the class 95% prediction interval, with named case studies",
       x = "Divergence time (My, log scale)",
       y = expression("Mitochondrial"~italic(K2P)~"(log scale)")) +
  theme_pub()
figS10     <- figS10_base
figS10_lab <- figS10_base +
  ggrepel::geom_text_repel(data = cases, aes(log_div, log_k2p, label = label),
    size = 3.0, lineheight = 0.9, fontface = "italic", colour = "grey10",
    box.padding = 0.9, min.segment.length = 0, seed = 1, max.overlaps = Inf)
save_fig(figS10,     "FigS7_mito_outliers",         width = 9.0, height = 6.4)
save_fig(figS10_lab, "FigS7_mito_outliers_labeled", width = 9.0, height = 6.4)

vp_m <- rd(file.path(MITO, "variance_partition.csv")) %>%
  mutate(Component = recode(component, pure_traits = "Pure traits", shared = "Shared",
                            pure_phylogeny = "Pure phylogeny", residual = "Residual"),
         Component = factor(Component, levels = names(COMP_PAL2)),
         ClassLab  = factor(unname(CLASS_LABELLER[class]), levels = c("Birds", "Mammals")),
         Model     = recode(model_type, avg_only = "Average", diff_only = "Difference"),
         Model     = factor(Model, levels = c("Average", "Difference")))
figS12 <- ggplot(vp_m, aes(ClassLab, median, fill = Component)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(median >= 0.04, sprintf("%.0f%%", 100 * median), "")),
            position = position_stack(vjust = 0.5), colour = "white", fontface = "bold", size = 3.1) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = COMP_PAL2, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.02))) +
  labs(title = "Variance partitioning of mitochondrial divergence",
       subtitle = "Median over 100 trees; average and difference models",
       x = NULL, y = "Share of variance") +
  theme_pub() + theme(panel.grid.major.x = element_blank())
save_fig(figS12, "FigS19_mito_variance_partition", width = 7.4, height = 5.2)

TRAIT_LABELS <- c(clutch = "Clutch size", mito_gc = "Mito GC content",
                  log_body_size = "Body size", log_gen_time = "Generation time",
                  log_pop_density = "Population density")
pretty_trait <- function(x) { o <- unname(TRAIT_LABELS[x])
  ifelse(is.na(o), str_to_sentence(str_replace_all(x, "_", " ")), o) }
MITO_BM <- c(Birds = MITO_PAL[["Aves"]], Mammals = MITO_PAL[["Mammalia"]])

load_w <- rd(file.path(MITO, "ppca_loadings_mean.csv")) %>%
  filter(axis %in% c("ppc1", "ppc2")) %>%
  pivot_wider(names_from = axis, values_from = loading) %>%
  mutate(ClassLab = factor(unname(CLASS_LABELLER[class]), levels = c("Birds", "Mammals")),
         trait_lab = pretty_trait(trait))
figS13A_base <- ggplot(load_w, aes(ppc1, ppc2)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.35) +
  geom_segment(aes(0, 0, xend = ppc1, yend = ppc2, colour = ClassLab),
               arrow = arrow(length = unit(6, "pt"), type = "closed"), linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~ ClassLab) +
  scale_colour_manual(values = MITO_BM, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.22)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  coord_equal(clip = "off") +
  labs(title = "A  Mitochondrial trait loadings on phylogenetic axes",
       x = "pPC1 (leading global axis)", y = "pPC2 (global)") +
  theme_pub() + theme(panel.spacing.x = unit(1.5, "cm"))
mito_trait_labels_layer <- ggrepel::geom_text_repel(aes(label = trait_lab), colour = "black", size = 2.9,
    box.padding = 1.1, force = 6, min.segment.length = 0, segment.size = 0.3,
    segment.colour = "grey65", max.iter = 20000, seed = 2, max.overlaps = Inf, show.legend = FALSE)
figS13A     <- figS13A_base  # unlabeled (arrows only)
figS13A_lab <- figS13A_base + mito_trait_labels_layer  # with individual trait labels

ax <- rd(file.path(MITO, "ppca_axis_effects.csv")) %>%
  filter(model_type == "ppca_avg") %>%
  mutate(ClassLab  = factor(unname(CLASS_LABELLER[class]), levels = c("Birds", "Mammals")),
         credible  = as.logical(credible),
         predictor = str_remove(predictor, " avg"),
         predictor = factor(predictor, levels = c("Sympatric", "Divergence Time",
                            "pPC3 (local)", "pPC2 (global)", "pPC1 (global)")))
dd <- dodge_forest(ax, "predictor", "ClassLab", width = 0.6)
figS13B <- ggplot(dd$df, aes(median, .y, colour = ClassLab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 2.5, linetype = "dotted", colour = "grey75", linewidth = 0.4) +
  geom_linerange(aes(xmin = q025, xmax = q975), linewidth = 0.7) +
  geom_point(aes(shape = credible), size = 2.6, fill = "white", stroke = 0.9) +
  scale_y_continuous(breaks = dd$breaks, labels = dd$labels) +
  scale_colour_manual(values = MITO_BM, name = NULL) +
  SHAPE_CRED +
  labs(title = "B  Phylogenetic-axis effects on mitochondrial divergence",
       subtitle = "Average model, pooled over 100 trees; pPC1-3 are trait axes, controls below the dotted line",
       x = "Standardized effect on log K2P", y = NULL) +
  theme_pub()
save_fig(figS13A,     "FigS8A_mito_ppca_biplot",         width = 8.0, height = 4.6)
save_fig(figS13A_lab, "FigS8A_mito_ppca_biplot_labeled", width = 8.0, height = 4.6)
save_fig(figS13B,     "FigS8B_mito_ppca_axis_effects",   width = 8.0, height = 4.4)
save_fig(figS13A_lab / figS13B + plot_layout(heights = c(1, 0.9)), "FigS8_mito_ppca_full",           width = 8.6, height = 9.2)
save_fig(figS13A     / figS13B + plot_layout(heights = c(1, 0.9)), "FigS8_mito_ppca_full_unlabeled", width = 8.6, height = 9.2)

cmp <- rd(file.path(MITO, "mito_k2p_compare.csv")) %>%
  filter(keep == 1, !is.na(k2p_whole), !is.na(k2p_coding)) %>%
  mutate(ClassLab = factor(unname(CLASS_LABELLER[taxonomic_group]), levels = names(LAB_MITO)))
rval <- cor(cmp$k2p_whole, cmp$k2p_coding)
figS14 <- ggplot(cmp, aes(k2p_whole, k2p_coding, colour = ClassLab)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 1.1, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = LAB_MITO, name = NULL, drop = FALSE) +
  coord_equal() +
  annotate("text", x = min(cmp$k2p_whole), y = max(cmp$k2p_coding), hjust = 0, vjust = 1, size = 3.4,
           label = sprintf("Pearson r = %.3f\nn = %d pairs", rval, nrow(cmp))) +
  labs(title = "Coding vs whole-mitogenome divergence",
       subtitle = "Coding K2P tracks whole-mitogenome K2P; dashed line is 1:1",
       x = expression("Whole-mitogenome"~italic(K2P)),
       y = expression("Coding (13 PCG)"~italic(K2P))) +
  theme_pub()
save_fig(figS14, "FigS9_whole_vs_coding_k2p", width = 5.8, height = 5.4)

sp_member <- rd(MEMBER) %>%
  dplyr::select(taxonomic_group, species1, species2, in_nuclear_yardstick, in_mito_yardstick) %>%
  pivot_longer(c(species1, species2), values_to = "species") %>%
  group_by(taxonomic_group, species) %>%
  summarise(nuc = as.integer(any(in_nuclear_yardstick == 1)),
            mit = as.integer(any(in_mito_yardstick == 1)), .groups = "drop") %>%
  mutate(Membership = case_when(nuc == 1 & mit == 1 ~ "Both",
                                nuc == 1 & mit == 0 ~ "Genome only",
                                nuc == 0 & mit == 1 ~ "Mito only"))
counts <- sp_member %>% count(taxonomic_group, Membership, name = "n") %>%
  mutate(ClassLab = factor(unname(CLASS_LABELLER[taxonomic_group]),
                           levels = c("Birds", "Mammals", "Reptiles", "Amphibians")),
         Membership = factor(Membership, levels = c("Genome only", "Both", "Mito only")))
totals <- counts %>% group_by(ClassLab) %>% summarise(tot = sum(n), .groups = "drop")
figS15 <- ggplot(counts, aes(ClassLab, n, fill = Membership)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n >= 6, n, "")), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.0) +
  geom_text(data = totals, aes(ClassLab, tot, label = tot), inherit.aes = FALSE,
            vjust = -0.4, size = 3.0, colour = "grey25") +
  scale_fill_manual(values = MEMBER_PAL, name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0.08))) +
  labs(title = "Species overlap between the nuclear and mitochondrial datasets",
       subtitle = "Congeneric-pair species assigned to each dataset, by class",
       x = NULL, y = "Number of species") +
  theme_pub() + theme(panel.grid.major.x = element_blank())
save_fig(figS15, "FigS10_dataset_overlap", width = 7.4, height = 4.8)

message("\nSupplementary mito figures written to: ", normalizePath(FIG_DIR))

qq_df <- mito %>% filter(Class_final %in% c("Aves","Mammalia","Reptilia"), !is.na(std_residual))
figS4 <- ggplot(qq_df, aes(sample = std_residual, colour = ClassLab)) +
  stat_qq(size = 0.9, alpha = 0.6) + stat_qq_line(colour = "grey40") +
  facet_wrap(~ ClassLab, scales = "free") +
  scale_colour_manual(values = LAB_MITO, guide = "none") +
  labs(title = "Mitochondrial yardstick residual Q–Q plots",
       x = "Theoretical quantile", y = "Standardized residual") +
  theme_pub()
save_fig(figS4, "FigS4_mito_qq", width = 9.0, height = 3.6)

PAIRS <- read_csv(file.path(MITO, "yardstick_pairs_all.csv"), show_col_types = FALSE)
mito_scatter <- function(cls_sci, model = c("avg","diff")) {
  model <- match.arg(model)
  cc <- c(body="Body.size..kg.", gen="Generation.time..years.",
          clutch="Clutch.size", dens="Population.Density", mitogc="Genome.GC.content")
  d <- PAIRS %>% filter(Class_final == cls_sci, timetree_div <= 15, k2p > 0) %>%
    mutate(log_k2p = log(k2p), z_div = as.numeric(scale(log(timetree_div))))
  build <- function(nm, col) {
    v1 <- suppressWarnings(as.numeric(d[[paste0(col, "_sp1")]]))
    v2 <- suppressWarnings(as.numeric(d[[paste0(col, "_sp2")]]))
    if (nm %in% c("mitogc","clutch")) { a <- rowMeans(cbind(v1, v2), na.rm = TRUE); df <- abs(v1 - v2) }
    else { a <- rowMeans(cbind(log(v1), log(v2)), na.rm = TRUE); df <- abs(log(v1) - log(v2)) }
    if (model == "avg") as.numeric(scale(a)) else as.numeric(scale(df))
  }
  long <- tibble(log_k2p = d$log_k2p,
                 `Divergence time` = d$z_div,
                 `Body size` = build("body", cc[["body"]]),
                 `Generation time` = build("gen", cc[["gen"]]),
                 `Clutch size` = build("clutch", cc[["clutch"]]),
                 `Mito GC` = build("mitogc", cc[["mitogc"]]))
  if (model == "avg") long$`Population density` <- build("dens", cc[["dens"]])
  long <- long %>% pivot_longer(-log_k2p, names_to = "Predictor", values_to = "z") %>% filter(!is.na(z))
  ggplot(long, aes(z, exp(log_k2p))) +
    geom_point(size = 0.8, alpha = 0.5, colour = MITO_PAL[[cls_sci]], stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.6) +
    facet_wrap(~ Predictor, scales = "free_x") + scale_y_log10() +
    labs(title = sprintf("%s %s trait model (mitochondrial)", cls_sci,
                         ifelse(model=="avg","average","difference")),
         x = "Standardized predictor", y = expression("Mitochondrial"~italic(K2P)~"(log)")) +
    theme_pub()
}
save_fig(mito_scatter("Mammalia","avg"),  "FigS5_mito_mammalia_avg_scatter",   width = 8.5, height = 5.4)
save_fig(mito_scatter("Aves","avg"),      "FigS6_mito_aves_avg_scatter",       width = 8.5, height = 5.4)
save_fig(mito_scatter("Mammalia","diff"), "FigS17_mito_mammalia_diff_scatter", width = 8.5, height = 5.0)
save_fig(mito_scatter("Aves","diff"),     "FigS18_mito_aves_diff_scatter",     width = 8.5, height = 5.0)

axd <- rd(file.path(MITO, "ppca_axis_effects.csv")) %>%
  filter(model_type == "ppca_diff") %>%
  mutate(ClassLab = factor(unname(CLASS_LABELLER[class]), levels = c("Birds","Mammals")),
         credible = as.logical(credible),
         predictor = str_remove(predictor, " diff"),
         predictor = factor(predictor, levels = c("Sympatric","Divergence Time",
                            "pPC3 (local)","pPC2 (global)","pPC1 (global)")))
ddd <- dodge_forest(axd, "predictor", "ClassLab", width = 0.6)
figS20 <- ggplot(ddd$df, aes(median, .y, colour = ClassLab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 2.5, linetype = "dotted", colour = "grey75", linewidth = 0.4) +
  geom_linerange(aes(xmin = q025, xmax = q975), linewidth = 0.7) +
  geom_point(aes(shape = credible), size = 2.6, fill = "white", stroke = 0.9) +
  scale_y_continuous(breaks = ddd$breaks, labels = ddd$labels) +
  scale_colour_manual(values = MITO_BM, name = NULL) +
  SHAPE_CRED +
  labs(title = "Mitochondrial pPCA axis effects (difference model)",
       x = "Standardized effect on log K2P", y = NULL) +
  theme_pub()
save_fig(figS20, "FigS20_mito_ppca_diff_effects", width = 8.0, height = 4.4)
