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

### candidate models for AIC
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

### class-specific models
for (cl in c("Aves", "Mammalia")) {
  subset <- complete_data %>% filter(Class == cl)
  if (nrow(subset) > 100) {
    m <- lmer(
      log_k2p ~ z_timetree_div + z_body_size + z_clutch +
        z_pop_density + z_genome_size + z_intersection +
        (1|Order) + (1|sp1_idx) + (1|sp2_idx),
      data = subset,
      REML = FALSE
    )
    print(summary(m)$coefficients)
    print(performance::r2_nakagawa(m))
  }
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
