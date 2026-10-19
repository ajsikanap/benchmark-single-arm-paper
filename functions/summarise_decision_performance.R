
rej      <- function(p, alpha) mean(p < alpha, na.rm = TRUE)
rej_mcse <- function(p, alpha) {
  x <- (p < alpha)[!is.na(p)]
  sqrt(mean(x) * (1 - mean(x)) / length(x))
}


summarise_decision_performance <- function(res_null,
                                           res_alt,
                                           alpha            = 0.10,
                                           scenario_label   = NULL,
                                           registry_n_label = NULL) {
  
  # Rejection rate and its Monte Carlo SE (binomial), over non-missing p-values
  rej <- function(p, alpha) mean(p < alpha, na.rm = TRUE)
  rej_mcse <- function(p, alpha) {
    x <- (p < alpha)[!is.na(p)]
    sqrt(mean(x) * (1 - mean(x)) / length(x))
  }
  
  typeI <- res_null |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n_null = dplyr::n(),
      typeI_one_sample_lr_registry_1s      = rej(one_sample_lr_p_one_sided_better_registry, alpha),
      typeI_one_sample_lr_registry_1s_mcse = rej_mcse(one_sample_lr_p_one_sided_better_registry, alpha),
      typeI_one_sample_lr_matched_1s       = rej(one_sample_lr_p_one_sided_better_matched, alpha),
      typeI_one_sample_lr_matched_1s_mcse  = rej_mcse(one_sample_lr_p_one_sided_better_matched, alpha),
      typeI_two_sample_lr_registry_1s      = rej(two_sample_lr_p_one_sided_naive, alpha),
      typeI_two_sample_lr_registry_1s_mcse = rej_mcse(two_sample_lr_p_one_sided_naive, alpha),
      typeI_two_sample_lr_matched_1s       = rej(two_sample_lr_p_one_sided_matched, alpha),
      typeI_two_sample_lr_matched_1s_mcse  = rej_mcse(two_sample_lr_p_one_sided_matched, alpha),
      .groups = "drop"
    )
  
  power <- res_alt |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      power_one_sample_lr_registry_1s      = rej(one_sample_lr_p_one_sided_better_registry, alpha),
      power_one_sample_lr_registry_1s_mcse = rej_mcse(one_sample_lr_p_one_sided_better_registry, alpha),
      power_one_sample_lr_matched_1s       = rej(one_sample_lr_p_one_sided_better_matched, alpha),
      power_one_sample_lr_matched_1s_mcse  = rej_mcse(one_sample_lr_p_one_sided_better_matched, alpha),
      power_two_sample_lr_registry_1s      = rej(two_sample_lr_p_one_sided_naive, alpha),
      power_two_sample_lr_registry_1s_mcse = rej_mcse(two_sample_lr_p_one_sided_naive, alpha),
      power_two_sample_lr_matched_1s       = rej(two_sample_lr_p_one_sided_matched, alpha),
      power_two_sample_lr_matched_1s_mcse  = rej_mcse(two_sample_lr_p_one_sided_matched, alpha),
      .groups = "drop"
    )
  
  benchmarks <- res_alt |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n_alt                = dplyr::n(),
      mean_trial_z1_est    = mean(trial_z1,      na.rm = TRUE),
      mean_trial_z0_est    = mean(trial_z0,      na.rm = TRUE),
      mean_matched_est     = mean(matched,       na.rm = TRUE),
      mean_registry_est    = mean(registry,      na.rm = TRUE),
      mean_rmst_trial_z1   = mean(rmst_trial_z1, na.rm = TRUE),
      mean_rmst_trial_z0   = mean(rmst_trial_z0, na.rm = TRUE),
      mean_rmst_matched    = mean(rmst_matched,  na.rm = TRUE),
      mean_rmst_registry   = mean(rmst_registry, na.rm = TRUE),
      mean_deltaS_matched  = mean(deltaS_obs,          na.rm = TRUE),
      mean_deltaS_registry = mean(deltaS_registry_obs, na.rm = TRUE),
      mean_deltaS_oracle   = mean(deltaS_oracle_obs,   na.rm = TRUE),
      .groups = "drop"
    )
  
  benchmarks |>
    dplyr::left_join(power, by = "method") |>
    dplyr::left_join(typeI, by = "method") |>
    dplyr::mutate(
      scenario   = scenario_label,
      registry_n = registry_n_label,
      .before    = method
    )
}