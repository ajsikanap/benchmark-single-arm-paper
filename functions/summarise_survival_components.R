summarise_survival_components <- function(res_boot, return_flextable = FALSE) {
  
  out <- res_boot |>
    dplyr::group_by(method,source) |>
    dplyr::summarise(
      
      ## ---- Trial treated (z1) ----
      mean_z1   = mean(trial_z1, na.rm = TRUE),
      sd_z1     = sd(trial_z1, na.rm = TRUE),
      var_z1    = var(trial_z1, na.rm = TRUE),
      q25_z1    = quantile(trial_z1, 0.25, na.rm = TRUE),
      median_z1 = quantile(trial_z1, 0.50, na.rm = TRUE),
      q75_z1    = quantile(trial_z1, 0.75, na.rm = TRUE),
      
      ## ---- Trial counterfactual (z0) ----
      mean_z0   = mean(trial_z0, na.rm = TRUE),
      sd_z0     = sd(trial_z0, na.rm = TRUE),
      var_z0    = var(trial_z0, na.rm = TRUE),
      q25_z0    = quantile(trial_z0, 0.25, na.rm = TRUE),
      median_z0 = quantile(trial_z0, 0.50, na.rm = TRUE),
      q75_z0    = quantile(trial_z0, 0.75, na.rm = TRUE),
      
      ## ---- Matched ----
      mean_matched   = mean(matched, na.rm = TRUE),
      sd_matched     = sd(matched, na.rm = TRUE),
      var_matched    = var(matched, na.rm = TRUE),
      q25_matched    = quantile(matched, 0.25, na.rm = TRUE),
      median_matched = quantile(matched, 0.50, na.rm = TRUE),
      q75_matched    = quantile(matched, 0.75, na.rm = TRUE),
      
      ## ---- Registry ----
      mean_registry   = mean(registry, na.rm = TRUE),
      sd_registry     = sd(registry, na.rm = TRUE),
      var_registry    = var(registry, na.rm = TRUE),
      q25_registry    = quantile(registry, 0.25, na.rm = TRUE),
      median_registry = quantile(registry, 0.50, na.rm = TRUE),
      q75_registry    = quantile(registry, 0.75, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  if (return_flextable) {
    flextable::flextable(out)
  } else {
    out
  }
}
