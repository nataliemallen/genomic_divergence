### setup
setwd("/Users/natal/Documents/Purdue/Genomic divergence/dataset/batch2/final_dataset")

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(performance)
library(MuMIn)
library(car)

### load and prep data
traits <- read_excel("Master_vertebrate_traits_12-08-25nma_FINAL.xlsx")
divergence <- read_csv("Complete_divergence_12-08-25.csv")

colnames(traits) <- make.names(colnames(traits))
colnames(divergence) <- make.names(colnames(divergence))

excluded_species <- traits %>%
  filter(Exclude. %in% c("yes", "Yes")) %>%
  pull(accession)

divergence_clean <- divergence %>%
  filter(!genome1 %in% excluded_species,
         !genome2 %in% excluded_species)

trait_cols <- c(
  "Body.size..kg.", "Generation.time..years.", "Clutch.size",
  "Population.Density", "Genome.size", "Genome.repeat.content",
  "Genome.GC.content", "area_km2"
)

traits <- traits %>%
  mutate(across(all_of(trait_cols), ~ suppressWarnings(as.numeric(.))))

traits_g1 <- traits %>%
  select(accession, Class, Order, all_of(trait_cols)) %>%
  rename(Class_g1 = Class, Order_g1 = Order) %>%
  rename_with(~ paste0(.x, "_g1"), all_of(trait_cols))

traits_g2 <- traits %>%
  select(accession, Class, Order, all_of(trait_cols)) %>%
  rename(Class_g2 = Class, Order_g2 = Order) %>%
  rename_with(~ paste0(.x, "_g2"), all_of(trait_cols))

merged_data <- divergence_clean %>%
  left_join(traits_g1, by = c("genome1" = "accession")) %>%
  left_join(traits_g2, by = c("genome2" = "accession")) %>%
  mutate(
    Class_final = coalesce(Class_g1, Class_g2),
    Order_final = coalesce(Order_g1, Order_g2)
  )

for (col in trait_cols) {
  g1 <- paste0(col, "_g1")
  g2 <- paste0(col, "_g2")
  avg <- paste0(col, "_avg")
  
  merged_data[[avg]] <- rowMeans(merged_data[, c(g1, g2)], na.rm = TRUE)
  merged_data[[avg]][is.nan(merged_data[[avg]])] <- NA
}

### analysis dataset
analysis_data <- merged_data %>%
  filter(Class_final %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"),
         !is.na(k2p),
         !is.na(timetree_div)) %>%
  mutate(
    log_k2p = log(k2p + 1e-8),
    log_timetree_div = log10(timetree_div + 1e-10),
    log_body_size = log10(Body.size..kg._avg + 0.001),
    log_genome_size = log10(Genome.size_avg + 1e-8),
    log_pop_density = log10(Population.Density_avg + 0.001),
    log_intersection = log10(intersection_km2 + 1),
    
    z_timetree_div = scale(log_timetree_div)[,1],
    z_body_size = scale(log_body_size)[,1],
    z_clutch = scale(Clutch.size_avg)[,1],
    z_pop_density = scale(log_pop_density)[,1],
    z_genome_size = scale(log_genome_size)[,1],
    z_intersection = scale(log_intersection)[,1],
    
    sp1_idx = factor(genome1),
    sp2_idx = factor(genome2),
    Class = factor(Class_final),
    Order = factor(Order_final)
  )

### candidate models for AIC with divergence time as fixed effect
complete_data <- analysis_data %>%
  filter(complete.cases(k2p, z_timetree_div, z_body_size, z_clutch,
                        z_pop_density, z_genome_size, z_intersection,
                        Order, sp1_idx, sp2_idx))

models <- list(
  m1 = lm(log_k2p ~ z_timetree_div, data = complete_data),
  
  m2 = lmer(log_k2p ~ z_timetree_div + (1|Order),
            data = complete_data, REML = FALSE),
  
  m3 = lmer(log_k2p ~ z_timetree_div +
              (1|Order) + (1|sp1_idx) + (1|sp2_idx),
            data = complete_data, REML = FALSE,
            control = lmerControl(optimizer = "bobyqa")),
  
  m4 = lmer(log_k2p ~ z_timetree_div + z_body_size + z_clutch +
              z_pop_density + z_genome_size + z_intersection +
              (1|Order) + (1|sp1_idx) + (1|sp2_idx),
            data = complete_data, REML = FALSE),
  
  m5 = lmer(log_k2p ~ z_timetree_div + z_body_size + z_clutch +
              z_pop_density + z_genome_size + z_intersection + Class +
              (1|Order) + (1|sp1_idx) + (1|sp2_idx),
            data = complete_data, REML = FALSE),
  
  m6 = lmer(log_k2p ~ z_timetree_div * Class + z_body_size + z_clutch +
              z_pop_density + z_genome_size + z_intersection +
              (1|Order) + (1|sp1_idx) + (1|sp2_idx),
            data = complete_data, REML = FALSE)
)

model_aic <- tibble(
  Model = names(models),
  AIC = map_dbl(models, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    weight = exp(-0.5*delta_AIC) / sum(exp(-0.5*delta_AIC))
  )

print(model_aic)

best_name <- model_aic$Model[1]
best_model <- models[[best_name]]

### model summary
print(summary(best_model))

r2_vals <- performance::r2_nakagawa(best_model)
print(r2_vals)

### GAM non-linearity test
if (nrow(complete_data) > 20) {
  gam_model <- gam(
    log_k2p ~ s(log_timetree_div) + s(log_body_size) +
      s(log_genome_size) + s(log_pop_density) +
      s(Order, bs = "re") + s(sp1_idx, bs = "re") + s(sp2_idx, bs = "re"),
    data = complete_data,
    method = "REML"
  )
  print(summary(gam_model))
}


### candidate models for AIC with k2p scaled by divergence time

# regress out divergence time
fit_divtime <- lm(log_k2p ~ z_timetree_div, data = complete_data)
complete_data$log_k2p_adj <- residuals(fit_divtime)

cat("Divergence time regressed out; residuals stored as log_k2p_adj\n")

# candidate models

models_list_adj <- list()

# intercept only (random effects only)
models_list_adj[["m1"]] <- lmer(log_k2p_adj ~ 1 +
                                  (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                data = complete_data, REML = FALSE)

# body size
models_list_adj[["m2"]] <- lmer(log_k2p_adj ~ z_body_size +
                                  (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                data = complete_data, REML = FALSE)

# body size + clutch
models_list_adj[["m3"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch +
                                  (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                data = complete_data, REML = FALSE)

# body size + clutch + population density
models_list_adj[["m4"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch + z_pop_density +
                                  (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                data = complete_data, REML = FALSE)

# body size + clutch + population density + genome size
models_list_adj[["m5"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch + z_pop_density + z_genome_size +
                                  (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                data = complete_data, REML = FALSE)

# all predictors including range overlap
models_list_adj[["m6"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch + z_pop_density + z_genome_size + z_intersection +
                                  (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                data = complete_data, REML = FALSE,
                                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

# compute AIC
model_aic_adj <- tibble(
  Model = names(models_list_adj),
  AIC = map_dbl(models_list_adj, ~ AIC(.x))
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  )

print(model_aic_adj)

# best model
best_model_name_adj <- model_aic_adj$Model[1]
best_model_adj <- models_list_adj[[best_model_name_adj]]

cat("\nBest model using divergence-adjusted k2p:", best_model_name_adj, "\n\n")
print(summary(best_model_adj))

### mammals and birds with k2p scaled by divergence time

for (cl in c("Aves", "Mammalia")) {
  cat("CLASS:", cl, "\n")
  
  class_data <- complete_data %>% filter(Class == cl)
  
  if (nrow(class_data) < 50) {
    next
  }
  
  # candidate models
  models_list_class <- list()
  
  models_list_class[["m1"]] <- lmer(log_k2p_adj ~ 1 +
                                      (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                    data = class_data, REML = FALSE)
  
  models_list_class[["m2"]] <- lmer(log_k2p_adj ~ z_body_size +
                                      (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                    data = class_data, REML = FALSE)
  
  models_list_class[["m3"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch +
                                      (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                    data = class_data, REML = FALSE)
  
  models_list_class[["m4"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch + z_pop_density +
                                      (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                    data = class_data, REML = FALSE)
  
  models_list_class[["m5"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch + z_pop_density + z_genome_size +
                                      (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                    data = class_data, REML = FALSE)
  
  models_list_class[["m6"]] <- lmer(log_k2p_adj ~ z_body_size + z_clutch + z_pop_density + z_genome_size + z_intersection +
                                      (1|Order) + (1|sp1_idx) + (1|sp2_idx),
                                    data = class_data, REML = FALSE,
                                    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
  # calculate AIC
  model_aic_class <- tibble(
    Model = names(models_list_class),
    AIC = map_dbl(models_list_class, ~ AIC(.x))
  ) %>%
    arrange(AIC) %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
    )
  
  print(model_aic_class)
  
  # best model
  best_model_name_class <- model_aic_class$Model[1]
  best_model_class <- models_list_class[[best_model_name_class]]
  
  cat("\nBest model for", cl, ":", best_model_name_class, "\n")
  print(summary(best_model_class))
}


### plots

# effect sizes
coef_df <- as.data.frame(summary(best_model)$coefficients) %>%
  rownames_to_column("Variable") %>%
  filter(Variable != "(Intercept)") %>%
  mutate(
    CI_lower = Estimate - 1.96*`Std. Error`,
    CI_upper = Estimate + 1.96*`Std. Error`
  )

p <- ggplot(coef_df, aes(Estimate, reorder(Variable, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  theme_minimal()

ggsave("figures/effect_sizes.png", p, width = 9, height = 6, dpi = 300)

p2 <- ggplot(analysis_data, aes(timetree_div, k2p, color = Class)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal()

ggsave("figures/divergence_vs_time.png", p2, width = 9, height = 6, dpi = 300)

##### scatterplots

dir.create("figures/scatterplots", showWarnings = FALSE)

plot_vars <- list(
  "Generation.time..years._avg" = "Generation time (years)",
  "Body.size..kg._avg"         = "Body size (kg)",
  "Clutch.size_avg"            = "Clutch size",
  "timetree_div"               = "Timetree divergence (MYA)",
  "Population.Density_avg"     = "Population density",
  "intersection_km2"           = "Range overlap (km²)"
)

for (v in names(plot_vars)) {
  
  label <- plot_vars[[v]]
  
  # apply log10 where appropriate
  if (all(analysis_data[[v]] > 0, na.rm = TRUE)) {
    analysis_data[[paste0(v, "_plot")]] <- log10(analysis_data[[v]])
    xvar <- paste0(v, "_plot")
    xlab <- paste0(label, " (log₁₀ scale)")
  } else {
    xvar <- v
    xlab <- label
  }
  
  p <- ggplot(analysis_data, aes_string(x = xvar, y = "k2p", color = "Class")) +
    geom_point(alpha = 0.4, size = 1.3) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_y_log10() +
    facet_wrap(~ Class, scales = "free") +
    theme_minimal(base_size = 13) +
    labs(
      x = xlab,
      y = "k2p divergence (log scale)",
      title = paste("k2p vs", label, "by Class")
    )
  
  outfile <- paste0("figures/scatterplots/k2p_vs_", v, ".png")
  ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

# plots of raw data

dir.create("figures/scatterplots_raw", showWarnings = FALSE)

raw_plot_vars <- list(
  "Generation.time..years._avg" = "Generation time (years)",
  "Body.size..kg._avg"         = "Body size (kg)",
  "Clutch.size_avg"            = "Clutch size",
  "timetree_div"               = "Timetree divergence (MYA)",
  "Population.Density_avg"     = "Population density",
  "intersection_km2"           = "Range overlap (km²)"
)

for (v in names(raw_plot_vars)) {
  
  label <- raw_plot_vars[[v]]
  
  p_raw <- ggplot(analysis_data, aes_string(x = v, y = "k2p", color = "Class")) +
    geom_point(alpha = 0.4, size = 1.3) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ Class, scales = "free") +
    theme_minimal(base_size = 13) +
    labs(
      x = label,
      y = "k2p divergence (raw scale)",
      title = paste("Raw-scale k2p vs", label, "by Class")
    )
  
  outfile_raw <- paste0("figures/scatterplots_raw/k2p_vs_raw_", v, ".png")
  ggsave(outfile_raw, p_raw, width = 10, height = 7, dpi = 300)
}
