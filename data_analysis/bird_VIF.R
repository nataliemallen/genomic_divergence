# first run data loading steps in bird DPMM

library(car)

compute_vif <- function(df, predictors) {
  
  missing <- setdiff(predictors, colnames(df))
  
  if (length(missing) > 0) {
    stop("Missing columns in dataframe: ", paste(missing, collapse = ", "))
  }
  
  df_sub <- df %>%
    select(all_of(predictors)) %>%
    drop_na()
  
  df_scaled <- as.data.frame(scale(df_sub))
  
  form <- as.formula(
    paste("y ~", paste(colnames(df_scaled), collapse = " + "))
  )
  
  df_scaled$y <- rnorm(nrow(df_scaled))
  
  model <- lm(form, data = df_scaled)
  
  car::vif(model)
}

bird_avg_predictors <- c(
  "z_timetree",
  "z_body_size_avg",
  "z_gen_time_avg",
  "z_clutch_avg",
  "z_pop_density_avg",
  "z_genome_size_avg",
  "z_gc_avg",
  "z_n50_avg"
)

bird_diff_predictors <- c(
  "z_timetree",
  "z_body_size_diff",
  "z_gen_time_diff",
  "z_clutch_diff",
  "z_genome_size_diff",
  "z_gc_diff",
  "z_n50_diff"
)

# df = your modeling dataframe BEFORE fitting

analysis_data_scaled <- scale_predictors(analysis_data)

vif_avg  <- compute_vif(analysis_data_scaled, bird_avg_predictors)
vif_diff <- compute_vif(analysis_data_scaled, bird_diff_predictors)

print("vif avg:")
print(vif_avg)
print("vif diff:")
print(vif_diff)
