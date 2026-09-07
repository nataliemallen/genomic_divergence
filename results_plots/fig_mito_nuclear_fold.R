# mito:nuclear divergence fold ratio per class suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
}))

if (file.exists("figures_setup.R")) source("figures_setup.R")

if (!exists("BASE"))  BASE  <- "."
if (!exists("YARD"))  YARD  <- file.path(BASE, "07_13_yardstick_results")
if (!exists("MITO"))  MITO  <- file.path(BASE, "mito_pipeline/results")
rd <- function(p) readr::read_csv(p, show_col_types = FALSE)
if (!exists("CLASS_PAL"))
  CLASS_PAL <- c(Aves = "#D55E00", Mammalia = "#0072B2", Reptilia = "#009E73")
if (!exists("theme_pub")) theme_pub <- function() theme_minimal(base_size = 13)
if (!exists("save_fig"))
  save_fig <- function(p, name, width = 7, height = 4.6) {
    dir.create("figures_extra", showWarnings = FALSE)
    for (ext in c("png", "pdf"))
      ggsave(file.path("figures_extra", paste0(name, ".", ext)),
             p, width = width, height = height, dpi = 300)
  }

ref_times <- c(1, 2, 5, 10, 15)

read_pred <- function(path, comp) {
  rd(path) %>%
    transmute(class, div_time_mya,
              expected_k2p = as.numeric(expected_k2p),
              extrapolation = as.logical(extrapolation),
              Compartment = comp)
}
nuc  <- read_pred(file.path(YARD, "reference_predictions_by_class.csv"), "Nuclear")
mito <- read_pred(file.path(MITO, "reference_predictions_by_class.csv"), "Mitochondrial")

fold <- inner_join(
  nuc  %>% select(class, div_time_mya, k2p_nuc = expected_k2p, extrap_nuc = extrapolation),
  mito %>% select(class, div_time_mya, k2p_mito = expected_k2p, extrap_mito = extrapolation),
  by = c("class", "div_time_mya")
) %>%
  filter(div_time_mya %in% ref_times,
         class %in% c("Aves", "Mammalia", "Reptilia"),
         !extrap_nuc, !extrap_mito) %>%
  mutate(fold = k2p_mito / k2p_nuc,
         class = factor(class, levels = c("Aves", "Mammalia", "Reptilia")),
         t = factor(div_time_mya, levels = ref_times, labels = paste0(ref_times, " My")))

fold_tbl <- fold %>%
  arrange(class, div_time_mya) %>%
  select(class, div_time_mya, k2p_nuc, k2p_mito, fold)
print(as.data.frame(fold_tbl), digits = 3)
write_csv(fold_tbl, "mito_nuclear_fold_by_class_time.csv")

p_fold <- ggplot(fold, aes(t, fold, colour = class, group = class)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.1f", fold)),
            vjust = -0.9, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = CLASS_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.08))) +
  labs(title = "Mitochondrial divergence relative to nuclear divergence",
       subtitle = "Fold difference (expected mito K2P ÷ expected nuclear K2P) by class",
       x = "Divergence time", y = "Mito : nuclear fold") +
  theme_pub()

save_fig(p_fold, "MitoNuclear_fold_by_class", width = 7.2, height = 4.8)
print(p_fold)

p_bar <- ggplot(fold, aes(t, fold, fill = class)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", fold)),
            position = position_dodge(0.75), vjust = -0.4, size = 3) +
  scale_fill_manual(values = CLASS_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.08))) +
  labs(title = "Mitochondrial : nuclear divergence fold by class and time",
       x = "Divergence time", y = "Mito : nuclear fold") +
  theme_pub()
save_fig(p_bar, "MitoNuclear_fold_by_class_bars", width = 7.2, height = 4.8)
print(p_bar)