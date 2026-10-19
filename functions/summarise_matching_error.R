summarise_matching_error <- function(res_boot,
                                     group_vars = c("method"),
                                     return_flextable = FALSE) {
  
  out <- res_boot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      
      n_iter = dplyr::n(),
      
      ## Survival matching error
      mean_err_S_matched_vs_trialcf = mean(err_S0_matched_vs_trialcf, na.rm = TRUE),
      
      mae_err_S_matched_vs_trialcf =  mean(abs(err_S0_matched_vs_trialcf), na.rm = TRUE),
      
      rmse_err_S_matched_vs_trialcf = sqrt(mean(err_S0_matched_vs_trialcf^2, na.rm = TRUE)),
      
      # Negativo per allineare la direzione: registry_true - trial_true_z0
      # coerente con matched_true_z0 - trial_true_z0
      mean_err_S_registry_vs_trialcf = mean(-err_S0_trial_vs_registry, na.rm = TRUE),
      
      mae_err_S_registry_vs_trialcf = mean(abs(err_S0_trial_vs_registry), na.rm = TRUE),
      
      rmse_err_S_registry_vs_trialcf =sqrt(mean(err_S0_trial_vs_registry^2, na.rm = TRUE)),
      
      ## RMST matching error
      mean_err_RMST_matched_vs_trialcf =mean(err_RMST0_matched_vs_trialcf, na.rm = TRUE),
      
      mae_err_RMST_matched_vs_trialcf = mean(abs(err_RMST0_matched_vs_trialcf), na.rm = TRUE),
      
      rmse_err_RMST_matched_vs_trialcf =sqrt(mean(err_RMST0_matched_vs_trialcf^2, na.rm = TRUE)),
      
      mean_err_RMST_registry_vs_trialcf = mean(-err_RMST0_trial_vs_registry, na.rm = TRUE),
      
      mae_err_RMST_registry_vs_trialcf =mean(abs(err_RMST0_trial_vs_registry), na.rm = TRUE),
      
      rmse_err_RMST_registry_vs_trialcf = sqrt(mean(err_RMST0_trial_vs_registry^2, na.rm = TRUE)),
      
      .groups = "drop"
    )
  
  if (return_flextable) flextable::flextable(out) else out
}
