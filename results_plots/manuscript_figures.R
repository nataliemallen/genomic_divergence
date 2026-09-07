# main-text figures 
setwd("")
source("figures_setup.R")

FIG_DIR <- file.path(BASE, "figures_main"); dir.create(FIG_DIR, showWarnings = FALSE)
rd <- function(p) readr::read_csv(p, show_col_types = FALSE)

COMP_PAL2 <- c(`Pure traits` = "#FF7D00", `Shared` = "#15616D",
               `Pure phylogeny` = "#E7DFC6", `Residual` = "#392C3A")

# axis breaks
XB <- log(c(0.1, 0.3, 1, 3, 10, 30)); XL <- c("0.1","0.3","1","3","10","30")
YBN <- log(c(0.001, 0.003, 0.01, 0.03, 0.1));       YLN <- c("0.001","0.003","0.01","0.03","0.1")
YBM <- log(c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3));  YLM <- c("0.001","0.003","0.01","0.03","0.1","0.3")

SHAPE_CRED <- scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                                 labels = c(`TRUE` = "credible", `FALSE` = "not credible"),
                                 name = NULL)

nuc  <- build_pairs()      %>% add_class_label()
mito <- build_mito_pairs() %>% add_class_label()

### 1A nuclear yardstick
ct_n   <- rd(file.path(YARD, "class_yardstick_table.csv"))
grid_n <- make_grid(ct_n, pal = LAB_NUC)
fig1A <- ggplot() +
  geom_ribbon(data = grid_n, aes(log_div, ymin = lwr, ymax = upr, fill = ClassLab), alpha = 0.12) +
  geom_point(data = nuc, aes(log_div, log_k2p, colour = ClassLab), size = 0.9, alpha = 0.35, stroke = 0) +
  geom_line(data = grid_n, aes(log_div, fit, colour = ClassLab), linewidth = 1.1) +
  scale_colour_manual(values = LAB_NUC, name = NULL, drop = FALSE) +
  scale_fill_manual(values = LAB_NUC, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = XB, labels = XL) +
  scale_y_continuous(breaks = YBN, labels = YLN) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.4))) +
  labs(title = "A  Nuclear genomic yardstick",
       x = "Divergence time (My, log scale)",
       y = expression("Genomic"~italic(K2P)~"(log scale)")) +
  theme_pub()

### 1B mitochondrial yardstick
ct_m   <- rd(file.path(MITO, "class_yardstick_table.csv"))
grid_m <- make_grid(ct_m, pal = LAB_MITO)
fig1B <- ggplot() +
  geom_ribbon(data = grid_m, aes(log_div, ymin = lwr, ymax = upr, fill = ClassLab), alpha = 0.12) +
  geom_point(data = mito, aes(log_div, log_k2p, colour = ClassLab), size = 0.9, alpha = 0.35, stroke = 0) +
  geom_line(data = grid_m, aes(log_div, fit, colour = ClassLab), linewidth = 1.1) +
  scale_colour_manual(values = LAB_MITO, name = NULL, drop = FALSE) +
  scale_fill_manual(values = LAB_MITO, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = XB, labels = XL) +
  scale_y_continuous(breaks = YBM, labels = YLM) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.4))) +
  labs(title = "B  Mitochondrial yardstick",
       x = "Divergence time (My, log scale)",
       y = expression("Mitochondrial"~italic(K2P)~"(log scale)")) +
  theme_pub()

### 1C reference predictions
ref_times <- c(1, 2, 5, 10, 15)
refC <- bind_rows(
  rd(file.path(YARD, "reference_predictions_by_class.csv")) %>% mutate(Compartment = "Nuclear"),
  rd(file.path(MITO, "reference_predictions_by_class.csv")) %>% mutate(Compartment = "Mitochondrial")) %>%

  filter(div_time_mya %in% ref_times, class %in% c("Aves", "Mammalia", "Reptilia"),
         !as.logical(extrapolation)) %>%
  mutate(ClassComp = factor(paste(class, Compartment, sep = "."), levels = names(CC_PAL)),
         ClassLab  = factor(unname(CLASS_LABELLER[class]), levels = c("Birds","Mammals","Reptiles")),
         Compartment = factor(Compartment, levels = c("Nuclear", "Mitochondrial")),
         t = factor(div_time_mya, levels = ref_times, labels = paste0(ref_times, " My")))
pdC <- position_dodge(width = 0.6)
fig1C <- ggplot(refC, aes(t, expected_k2p, colour = ClassComp, group = ClassComp)) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = pdC, linewidth = 0.8) +
  geom_point(position = pdC, size = 2.2) +
  facet_grid(Compartment ~ ., switch = "y") +
  scale_colour_manual(values = CC_PAL, name = NULL, labels = function(x) gsub("\\.", " ", x)) +
  scale_y_log10() +
  labs(title = "C  Expected divergence at reference times",
       subtitle = "Predicted K2P (95% prediction interval)",
       x = NULL, y = expression("Expected"~italic(K2P)~"(log scale)")) +
  theme_pub() +
  theme(legend.text = element_text(size = rel(0.75)))

### 1D residuals by Order
ord_res <- bind_rows(
  nuc  %>% mutate(Compartment = "Nuclear"),
  mito %>% mutate(Compartment = "Mitochondrial")) %>%
  filter(Class_final %in% c("Aves", "Mammalia", "Reptilia"), !is.na(Order_final)) %>%
  add_count(Compartment, Class_final, Order_final, name = "n_ord") %>%
  filter(n_ord >= 5) %>%
  mutate(ClassComp   = factor(paste(Class_final, Compartment, sep = "."), levels = names(CC_PAL)),
         Compartment = factor(Compartment, levels = c("Nuclear", "Mitochondrial")),
         Order_final = fct_reorder(Order_final, residual, .fun = median))
fig1D <- ggplot(ord_res, aes(residual, Order_final, colour = ClassComp, fill = ClassComp)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_boxplot(alpha = 0.20, outlier.size = 0.5, linewidth = 0.4,
               position = position_dodge2(preserve = "single")) +
  scale_colour_manual(values = CC_PAL, name = NULL, labels = function(x) gsub("\\.", " ", x)) +
  scale_fill_manual(values = CC_PAL, guide = "none") +
  facet_grid(ClassLab ~ ., scales = "free_y", space = "free_y") +
  labs(title = "D  Residual divergence by Order",
       subtitle = "Deviation from the class expectation (log K2P units); 0 = on the yardstick",
       x = "Residual (faster → / ← slower than expected)", y = NULL) +
  theme_pub() +
  theme(legend.text = element_text(size = rel(0.75)), panel.grid.major.y = element_blank())

save_fig(fig1A, "Fig1A_nuclear_yardstick",   width = 5.6, height = 4.6)
save_fig(fig1B, "Fig1B_mito_yardstick",       width = 5.6, height = 4.6)
save_fig(fig1C, "Fig1C_reference_times",      width = 7.0, height = 4.6)
save_fig(fig1D, "Fig1D_residuals_by_order",   width = 7.4, height = 7.6)
fig1 <- (fig1A | fig1B) / (fig1C | fig1D) + plot_layout(heights = c(1, 1.15))
save_fig(fig1, "Fig1_yardstick_full", width = 13.5, height = 12.0)

### 2 nuclear yardstick outliers
flagged <- nuc %>%
  mutate(flag = case_when(std_residual >  1.96 ~ "Higher than expected",
                          std_residual < -1.96 ~ "Lower than expected",
                          TRUE ~ "On the yardstick"))
pick <- function(g1, g2) nuc %>% filter((grepl(g1, sp1) & grepl(g2, sp2)) |
                                        (grepl(g2, sp1) & grepl(g1, sp2)))
cases <- bind_rows(
  pick("pteronotus_mesoamericanus", "pteronotus_parnellii") %>% mutate(label = "Pteronotus mesoamericanus / P. parnellii"),
  pick("gavia_immer", "gavia_adamsii")                       %>% mutate(label = "Gavia immer / G. adamsii"),
  pick("pyrrhura_roseifrons", "pyrrhura_egregia")            %>% mutate(label = "Pyrrhura roseifrons / P. egregia"),
  pick("alouatta_juara", "alouatta_macconnelli")             %>% mutate(label = "Alouatta juara / A. macconnelli")
) %>% distinct()
fig2_base <- ggplot() +
  geom_line(data = grid_n, aes(log_div, fit, group = ClassLab), colour = "grey70", linewidth = 0.7) +
  geom_point(data = filter(flagged, flag == "On the yardstick"),
             aes(log_div, log_k2p), colour = "grey80", size = 0.8, alpha = 0.5, stroke = 0) +
  geom_point(data = filter(flagged, flag != "On the yardstick"),
             aes(log_div, log_k2p, colour = flag), size = 2.0, alpha = 0.9, stroke = 0) +
  scale_colour_manual(values = OUT_COL, name = NULL) +
  geom_point(data = cases, aes(log_div, log_k2p), shape = 21, fill = CASE_FILL,
             colour = "grey15", size = 3.4, stroke = 0.7) +
  scale_x_continuous(breaks = XB, labels = XL) +
  scale_y_continuous(breaks = YBN, labels = YLN) +
  labs(title = "The yardstick flags biologically meaningful outliers",
       subtitle = "Pairs outside the class 95% prediction interval, with named case studies",
       x = "Divergence time (My, log scale)",
       y = expression("Genomic"~italic(K2P)~"(log scale)")) +
  theme_pub()
# without and with case-study labels 
fig2 <- fig2_base
fig2_lab <- fig2_base +
  ggrepel::geom_text_repel(data = cases, aes(log_div, log_k2p, label = label),
    size = 3.0, lineheight = 0.9, fontface = "italic", colour = "grey10",
    box.padding = 0.9, min.segment.length = 0, seed = 1, max.overlaps = Inf)
save_fig(fig2,     "Fig2_outliers",           width = 9.0, height = 6.4)
save_fig(fig2_lab, "Fig2_outliers_labeled",   width = 9.0, height = 6.4)

### 3 PMM results
PRED_ORDER <- c("Divergence Time", "Body Size", "Generation Time", "Clutch Size",
                "Population Density", "Sympatric",  # shared
                "GC Content", "Genome Size", "N50",  # nuclear-only
                "Mito GC")  # mitochondrial-only
pe <- bind_rows(
  rd(file.path(PMM,  "pooled_effects_all_classes_models.csv")) %>% mutate(Compartment = "Nuclear"),
  rd(file.path(MITO, "pooled_effects_all_classes_models.csv")) %>% mutate(Compartment = "Mitochondrial")) %>%
  filter(model_type == "avg_only", class %in% c("birds", "mammals"), predictor %in% PRED_ORDER) %>%
  mutate(ClassSci    = ifelse(class == "birds", "Aves", "Mammalia"),
         ClassLab    = factor(unname(CLASS_LABELLER[class]), levels = c("Birds", "Mammals")),
         Compartment = factor(Compartment, levels = c("Nuclear", "Mitochondrial")),
         ClassComp   = factor(paste(ClassSci, Compartment, sep = "."), levels = names(CC_PAL)),
         predictor   = factor(predictor, levels = rev(PRED_ORDER)),
         credible    = as.logical(credible))

dd3 <- dodge_forest(pe, "predictor", "Compartment", width = 0.6)
fig3 <- ggplot(dd3$df, aes(median, .y, colour = ClassComp)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_linerange(aes(xmin = q025, xmax = q975), linewidth = 0.7) +
  geom_point(aes(shape = credible), size = 2.5, fill = "white", stroke = 0.9) +
  facet_wrap(~ ClassLab) +
  scale_y_continuous(breaks = dd3$breaks, labels = dd3$labels) +
  scale_colour_manual(values = CC_PAL, name = NULL, labels = function(x) gsub("\\.", " ", x)) +
  SHAPE_CRED +
  labs(title = "Trait effects on genomic divergence: nuclear vs. mitochondrial",
       subtitle = "Average-model standardized coefficients (95% CrI, pooled over 100 trees); lower rows are compartment-specific predictors",
       x = "Standardized effect on log K2P", y = NULL) +
  theme_pub()
save_fig(fig3, "Fig3_pmm_coefficients", width = 10.0, height = 6.0)

### 4 variance partitioning
names(COMP_PAL2) <- c(
  "Pure trait variance",
  "Phylogenetically structured trait variance",
  "Pure phylogeny variance",
  "Residual variance"
)

vp <- bind_rows(
  rd(file.path(PMM,  "variance_partition.csv")) %>% mutate(Compartment = "Nuclear"),
  rd(file.path(MITO, "variance_partition.csv")) %>% mutate(Compartment = "Mitochondrial")
) %>%
  filter(model_type == "avg_only", class %in% c("birds", "mammals")) %>%
  mutate(
    Component = recode(
      component,
      pure_traits    = "Pure trait variance",
      shared         = "Phylogenetically structured trait variance",
      pure_phylogeny = "Pure phylogeny variance",
      residual       = "Residual variance"
    ),
    Component   = factor(Component, levels = names(COMP_PAL2)),
    ClassLab    = factor(unname(CLASS_LABELLER[class]), levels = c("Birds", "Mammals")),
    Compartment = factor(Compartment, levels = c("Nuclear", "Mitochondrial"))
  )

fig4 <- ggplot(vp, aes(Compartment, median, fill = Component)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(median >= 0.04, sprintf("%.0f%%", 100 * median), "")),
            position = position_stack(vjust = 0.5), colour = "white", fontface = "bold", size = 3.0) +
  facet_wrap(~ ClassLab) +
  scale_fill_manual(values = COMP_PAL2, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.02))) +
  labs(title = "Traits act more independently of ancestry in mtDNA than in the nuclear genome",
       subtitle = "Partition of explained variance in divergence (average model, median over 100 trees)",
       x = NULL, y = "Share of variance") +
  theme_pub() + theme(panel.grid.major.x = element_blank())
save_fig(fig4, "Fig4_variance_partition", width = 9.0, height = 5.2)

### 5 ppca
TRAIT_LABELS <- c(clutch = "Clutch size", gc = "GC content", log_body_size = "Body size",
                  log_gen_time = "Generation time", log_genome_size = "Genome size",
                  log_n50 = "N50", log_pop_density = "Population density")
pretty_trait <- function(x) { o <- unname(TRAIT_LABELS[x])
  ifelse(is.na(o), str_to_sentence(str_replace_all(x, "_", " ")), o) }

### 5A loadings biplot
load_w <- rd(file.path(PMM, "ppca_loadings_mean.csv")) %>%
  filter(axis %in% c("ppc1", "ppc2")) %>%
  pivot_wider(names_from = axis, values_from = loading) %>%
  mutate(ClassLab = factor(unname(CLASS_LABELLER[class]), levels = c("Birds", "Mammals")),
         trait_lab = pretty_trait(trait))
fig5A_base <- ggplot(load_w, aes(ppc1, ppc2)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.35) +
  geom_segment(aes(0, 0, xend = ppc1, yend = ppc2, colour = ClassLab),
               arrow = arrow(length = unit(6, "pt"), type = "closed"), linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~ ClassLab) +
  scale_colour_manual(values = c(Birds = NUC_PAL[["Aves"]], Mammals = NUC_PAL[["Mammalia"]]), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.22)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  coord_equal(clip = "off") +
  labs(title = "A  Trait loadings on phylogenetic axes",
       subtitle = "pPCA global axes (Abouheif proximity); arrows are trait contributions",
       x = "pPC1 (leading global axis)", y = "pPC2 (global)") +
  theme_pub() + theme(panel.spacing.x = unit(1.5, "cm"))
trait_labels_layer <- ggrepel::geom_text_repel(aes(label = trait_lab), colour = "black", size = 2.9,
    box.padding = 1.1, force = 6, min.segment.length = 0, segment.size = 0.3,
    segment.colour = "grey65", max.iter = 20000, seed = 2, max.overlaps = Inf, show.legend = FALSE)
fig5A     <- fig5A_base  # unlabeled (arrows only)
fig5A_lab <- fig5A_base + trait_labels_layer  # with individual trait labels

### 5B axis effects 
axn <- rd(file.path(PMM,  "ppca_axis_effects.csv")) %>% filter(model_type == "ppca_avg") %>%
  transmute(class, predictor, est = estimate, lo = ci_lo, hi = ci_hi,
            credible = as.logical(significant), Compartment = "Nuclear")
axm <- rd(file.path(MITO, "ppca_axis_effects.csv")) %>% filter(model_type == "ppca_avg") %>%
  transmute(class, predictor, est = median, lo = q025, hi = q975,
            credible = as.logical(credible), Compartment = "Mitochondrial")
ax <- bind_rows(axn, axm) %>%
  filter(class %in% c("birds", "mammals")) %>%
  mutate(ClassSci    = ifelse(class == "birds", "Aves", "Mammalia"),
         Compartment = factor(Compartment, levels = c("Nuclear", "Mitochondrial")),
         ClassComp   = factor(paste(ClassSci, Compartment, sep = "."), levels = names(CC_PAL)),
         predictor   = str_remove(predictor, " avg"),
         predictor   = factor(predictor, levels = c("Sympatric", "Divergence Time",
                              "pPC3 (local)", "pPC2 (global)", "pPC1 (global)")))
dd5 <- dodge_forest(ax, "predictor", "ClassComp", width = 0.78)
fig5B <- ggplot(dd5$df, aes(est, .y, colour = ClassComp)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 2.5, linetype = "dotted", colour = "grey75", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6) +
  geom_point(aes(shape = credible), size = 2.3, fill = "white", stroke = 0.8) +
  scale_y_continuous(breaks = dd5$breaks, labels = dd5$labels) +
  scale_colour_manual(values = CC_PAL, name = NULL, labels = function(x) gsub("\\.", " ", x)) +
  SHAPE_CRED +
  labs(title = "B  Relationship of phylogenetic axes to divergence",
       subtitle = "Pooled standardized effects (average model); nuclear and mitochondrial; pPC1-3 are trait axes, controls below the dotted line",
       x = "Standardized effect on log K2P", y = NULL) +
  theme_pub()
save_fig(fig5A,     "Fig5A_ppca_biplot",           width = 8.0, height = 4.6)
save_fig(fig5A_lab, "Fig5A_ppca_biplot_labeled",   width = 8.0, height = 4.6)
save_fig(fig5B,     "Fig5B_ppca_axis_effects",     width = 8.0, height = 4.8)
save_fig(fig5A_lab / fig5B + plot_layout(heights = c(1, 0.9)), "Fig5_ppca_full",           width = 8.6, height = 9.4)
save_fig(fig5A     / fig5B + plot_layout(heights = c(1, 0.9)), "Fig5_ppca_full_unlabeled", width = 8.6, height = 9.4)
