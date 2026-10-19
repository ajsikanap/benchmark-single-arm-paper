summarise_boot_results_old <- function(res_boot, return_flextable = FALSE) {
  
  out <- res_boot |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      
      
      ## Cox model 0: source only
      beta_mean_0   = mean(beta_hat_0, na.rm = TRUE),
      bias_0        = mean(bias_beta_0, na.rm = TRUE),
      rmse_0        = sqrt(mean(bias_beta_0^2, na.rm = TRUE)),
      var_emp_0     = var(beta_hat_0, na.rm = TRUE),
      power_2sided_0 = mean(p_val_0, na.rm = TRUE),
      power_1sided_0 = mean(sign_pval_0, na.rm = TRUE),
      
      
      ## Cox model 1: source + prog_score
      
      beta_mean_1   = mean(beta_hat_1, na.rm = TRUE),
      bias_1        = mean(bias_beta_1, na.rm = TRUE),
      rmse_1        = sqrt(mean(bias_beta_1^2, na.rm = TRUE)),
      var_emp_1     = var(beta_hat_1, na.rm = TRUE),
      power_2sided_1 = mean(p_val_1, na.rm = TRUE),
      power_1sided_1 = mean(sign_pval_1, na.rm = TRUE),
      
      
      ## Cox model 2: source + score
      
      beta_mean_2   = mean(beta_hat_2, na.rm = TRUE),
      bias_2        = mean(bias_beta_2, na.rm = TRUE),
      rmse_2        = sqrt(mean(bias_beta_2^2, na.rm = TRUE)),
      var_emp_2     = var(beta_hat_2, na.rm = TRUE),
      power_2sided_2 = mean(p_val_2, na.rm = TRUE),
      power_1sided_2 = mean(sign_pval_2, na.rm = TRUE),
      
      ## Survival at tau = 12
      
      surv_12mean_trial_z1   = mean(trial_z1, na.rm = TRUE),
      surv_12mean_trial_z0   = mean(trial_z0, na.rm = TRUE),
      surv_12mean_matched    = mean(matched, na.rm = TRUE),
      surv_12mean_registry   = mean(registry, na.rm = TRUE),
      
      bias_S1_trial          = mean(bias_S1_trial, na.rm = TRUE),
      bias_S0_matched        = mean(bias_S0_matched, na.rm = TRUE),
      bias_S0_registry       = mean(bias_S0_registry, na.rm = TRUE),
      bias_deltaS_matched    = mean(bias_deltaS_matched, na.rm = TRUE),
      bias_deltaS_registry   = mean(bias_deltaS_registry, na.rm = TRUE),
      err_S0_matched_vs_trialcf = mean(err_S0_matched_vs_trialcf, na.rm=TRUE),
      err_S0_trial_vs_registry = mean(err_S0_trial_vs_registry, na.rm=TRUE),
      
      
      ## RMST up to tau = 12
      
      rmst_mean_trial_z1     = mean(rmst_trial_z1, na.rm = TRUE),
      rmst_mean_trial_z0     = mean(rmst_trial_z0, na.rm = TRUE),
      rmst_mean_matched      = mean(rmst_matched, na.rm = TRUE),
      rmst_mean_registry     = mean(rmst_registry, na.rm = TRUE),
      
      bias_RMST1_trial       = mean(bias_RMST1_trial, na.rm = TRUE),
      bias_RMST0_matched     = mean(bias_RMST0_matched, na.rm = TRUE),
      bias_RMST0_registry    = mean(bias_RMST0_registry, na.rm = TRUE),
      bias_deltaRMST_matched = mean(bias_deltaRMST_matched, na.rm = TRUE),
      bias_deltaRMST_registry = mean(bias_deltaRMST_registry, na.rm = TRUE),
      err_RMST0_matched_vs_trialcf = mean(err_RMST0_matched_vs_trialcf, na.rm=TRUE),
      err_RMST0_trial_vs_registry = mean(err_RMST0_trial_vs_registry, na.rm=TRUE),
      .groups = "drop"
    )
  
  if (return_flextable) {
    return(flextable::flextable(out))
  } else {
    return(out)
  }
}


summarise_boot_results <- function(res_boot, return_flextable = FALSE) {
  
  out <- res_boot |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      
      # ── Cox model 0: source only ──────────────────────────────────────────
      beta_mean_0    = mean(beta_hat_0,   na.rm = TRUE),
      beta_sd_0      = sd(beta_hat_0,     na.rm = TRUE),
      bias_0         = mean(bias_beta_0,  na.rm = TRUE),
      rmse_0         = sqrt(mean(bias_beta_0^2, na.rm = TRUE)),
      var_emp_0      = var(beta_hat_0,    na.rm = TRUE),
      power_2sided_0 = mean(p_val_0   < 0.05, na.rm = TRUE),
      power_1sided_0 = mean(sign_pval_0,       na.rm = TRUE),
      
      # ── Cox model 1: source + prog_score ─────────────────────────────────
      beta_mean_1    = mean(beta_hat_1,   na.rm = TRUE),
      beta_sd_1      = sd(beta_hat_1,     na.rm = TRUE),
      bias_1         = mean(bias_beta_1,  na.rm = TRUE),
      rmse_1         = sqrt(mean(bias_beta_1^2, na.rm = TRUE)),
      var_emp_1      = var(beta_hat_1,    na.rm = TRUE),
      power_2sided_1 = mean(p_val_1   < 0.05, na.rm = TRUE),
      power_1sided_1 = mean(sign_pval_1,       na.rm = TRUE),
      
      # ── Cox model 2: source + score ───────────────────────────────────────
      beta_mean_2    = mean(beta_hat_2,   na.rm = TRUE),
      beta_sd_2      = sd(beta_hat_2,     na.rm = TRUE),
      bias_2         = mean(bias_beta_2,  na.rm = TRUE),
      rmse_2         = sqrt(mean(bias_beta_2^2, na.rm = TRUE)),
      var_emp_2      = var(beta_hat_2,    na.rm = TRUE),
      power_2sided_2 = mean(p_val_2   < 0.05, na.rm = TRUE),
      power_1sided_2 = mean(sign_pval_2,       na.rm = TRUE),
      
      # ── Operating characteristics ─────────────────────────────────────────
      reject_one_sample_registry = mean(one_sample_lr_p_one_sided_better_registry < 0.10, na.rm = TRUE),
      reject_one_sample_matched  = mean(one_sample_lr_p_one_sided_better_matched  < 0.10, na.rm = TRUE),
      reject_two_sample_naive    = mean(two_sample_lr_p_one_sided_naive           < 0.10, na.rm = TRUE),
      reject_two_sample_matched  = mean(two_sample_lr_p_one_sided_matched         < 0.10, na.rm = TRUE),
      
      # ── Survival at tau ───────────────────────────────────────────────────
      surv_mean_trial_z1  = mean(trial_z1, na.rm = TRUE),
      surv_sd_trial_z1    = sd(trial_z1,   na.rm = TRUE),
      surv_mean_trial_z0  = mean(trial_z0, na.rm = TRUE),
      surv_sd_trial_z0    = sd(trial_z0,   na.rm = TRUE),
      surv_mean_matched   = mean(matched,  na.rm = TRUE),
      surv_sd_matched     = sd(matched,    na.rm = TRUE),
      surv_mean_registry  = mean(registry, na.rm = TRUE),
      surv_sd_registry    = sd(registry,   na.rm = TRUE),
      
      # Latent survival means
      surv_mean_trial_true_z1  = mean(trial_true_z1,  na.rm = TRUE),
      surv_mean_trial_true_z0  = mean(trial_true_z0,  na.rm = TRUE),
      surv_mean_matched_true_z0= mean(matched_true_z0,na.rm = TRUE),
      surv_mean_registry_true  = mean(registry_true,  na.rm = TRUE),
      
      # Survival biases
      bias_S1_trial        = mean(bias_S1_trial,        na.rm = TRUE),
      bias_S0_matched      = mean(bias_S0_matched,      na.rm = TRUE),
      bias_S0_registry     = mean(bias_S0_registry,     na.rm = TRUE),
      bias_deltaS_matched  = mean(bias_deltaS_matched,  na.rm = TRUE),
      bias_deltaS_registry = mean(bias_deltaS_registry, na.rm = TRUE),
      
      # Survival bias decomposition
      mean_err_trial_z1_est          = mean(err_trial_z1_est,          na.rm = TRUE),
      mean_err_matched_z0_est        = mean(err_matched_z0_est,        na.rm = TRUE),
      mean_err_registry_z0_est       = mean(err_registry_z0_est,       na.rm = TRUE),
      mean_err_matching_z0           = mean(err_matching_z0,           na.rm = TRUE),
      mean_err_registry_selection_z0 = mean(err_registry_selection_z0, na.rm = TRUE),
      mean_sampling_error_deltaS     = mean(sampling_error_deltaS,     na.rm = TRUE),
      mean_bias_deltaS_decomp_matched= mean(bias_deltaS_decomp_matched,na.rm = TRUE),
      mean_bias_deltaS_decomp_naive  = mean(bias_deltaS_decomp_naive,  na.rm = TRUE),
      
      # Matching error survival
      mean_err_S0_matched_vs_trialcf = mean(err_S0_matched_vs_trialcf, na.rm = TRUE),
      mean_err_S0_trial_vs_registry  = mean(err_S0_trial_vs_registry,  na.rm = TRUE),
      
      # ── RMST at tau ───────────────────────────────────────────────────────
      rmst_mean_trial_z1  = mean(rmst_trial_z1, na.rm = TRUE),
      rmst_sd_trial_z1    = sd(rmst_trial_z1,   na.rm = TRUE),
      rmst_mean_trial_z0  = mean(rmst_trial_z0, na.rm = TRUE),
      rmst_sd_trial_z0    = sd(rmst_trial_z0,   na.rm = TRUE),
      rmst_mean_matched   = mean(rmst_matched,  na.rm = TRUE),
      rmst_sd_matched     = sd(rmst_matched,    na.rm = TRUE),
      rmst_mean_registry  = mean(rmst_registry, na.rm = TRUE),
      rmst_sd_registry    = sd(rmst_registry,   na.rm = TRUE),
      
      # Latent RMST means
      rmst_mean_trial_true_z1   = mean(rmst_trial_true_z1,   na.rm = TRUE),
      rmst_mean_trial_true_z0   = mean(rmst_trial_true_z0,   na.rm = TRUE),
      rmst_mean_matched_true_z0 = mean(rmst_matched_true_z0, na.rm = TRUE),
      rmst_mean_registry_true   = mean(rmst_registry_true,   na.rm = TRUE),
      
      # RMST biases
      bias_RMST1_trial         = mean(bias_RMST1_trial,         na.rm = TRUE),
      bias_RMST0_matched       = mean(bias_RMST0_matched,       na.rm = TRUE),
      bias_RMST0_registry      = mean(bias_RMST0_registry,      na.rm = TRUE),
      bias_deltaRMST_matched   = mean(bias_deltaRMST_matched,   na.rm = TRUE),
      bias_deltaRMST_registry  = mean(bias_deltaRMST_registry,  na.rm = TRUE),
      
      # RMST bias decomposition
      mean_err_trial_z1_rmst_est        = mean(err_trial_z1_rmst_est,        na.rm = TRUE),
      mean_err_matched_z0_rmst_est      = mean(err_matched_z0_rmst_est,      na.rm = TRUE),
      mean_err_registry_z0_rmst_est     = mean(err_registry_z0_rmst_est,     na.rm = TRUE),
      mean_err_matching_rmst_z0         = mean(err_matching_rmst_z0,         na.rm = TRUE),
      mean_err_registry_selection_rmst  = mean(err_registry_selection_rmst,  na.rm = TRUE),
      mean_sampling_error_deltaRMST     = mean(sampling_error_deltaRMST,     na.rm = TRUE),
      mean_bias_deltaRMST_decomp_matched= mean(bias_deltaRMST_decomp_matched,na.rm = TRUE),
      mean_bias_deltaRMST_decomp_naive  = mean(bias_deltaRMST_decomp_naive,  na.rm = TRUE),
      
      # Matching error RMST
      mean_err_RMST0_matched_vs_trialcf = mean(err_RMST0_matched_vs_trialcf, na.rm = TRUE),
      mean_err_RMST0_trial_vs_registry  = mean(err_RMST0_trial_vs_registry,  na.rm = TRUE),
      
      # ── Balance diagnostics ───────────────────────────────────────────────
      score_smd_before         = mean(score_smd_before,         na.rm = TRUE),
      score_smd_matched        = mean(score_smd_matched,        na.rm = TRUE),
      score_smd_reduction_abs  = mean(score_smd_reduction_abs,  na.rm = TRUE),
      score_smd_reduction_rel  = mean(score_smd_reduction_rel,  na.rm = TRUE),
      score_var_ratio_before   = mean(score_var_ratio_before,   na.rm = TRUE),
      score_var_ratio_matched  = mean(score_var_ratio_matched,  na.rm = TRUE),
      score_pair_distance      = mean(score_pair_distance,      na.rm = TRUE),
      n_control_all            = mean(n_control_all,            na.rm = TRUE),
      n_trial_all              = mean(n_trial_all,              na.rm = TRUE),
      n_control_matched        = mean(n_control_matched,        na.rm = TRUE),
      n_control_ess            = mean(n_control_ess,            na.rm = TRUE),
      n_control_discarded      = mean(n_control_discarded,      na.rm = TRUE),
      prop_control_discarded   = mean(prop_control_discarded,   na.rm = TRUE),
      smd_cov_mean_matched     = mean(smd_cov_mean_matched,     na.rm = TRUE),
      smd_cov_max_matched      = mean(smd_cov_max_matched,      na.rm = TRUE),
      smd_cov_min_matched      = mean(smd_cov_min_matched,      na.rm = TRUE),
      smd_cov_prop_gt_0.1      = mean(smd_cov_prop_gt_0.1,      na.rm = TRUE),
      smd_cov_prop_gt_0.2      = mean(smd_cov_prop_gt_0.2,      na.rm = TRUE),
      
      .groups = "drop"
    )
  
  if (return_flextable) {
    return(flextable::flextable(out))
  } else {
    return(out)
  }
}