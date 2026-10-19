# make_na_row.R
# This function builds a row with all columns as NA, to be used when the
# matching fails. It does not include effect_setting and the p-value.

make_na_row <- function(iter, seed, method, tau,
                        beta_unadj_star, beta_adj_prog_star, beta_adj_prop_star,
                        S0_true_tau, S1_true_tau, deltaS_true,
                        RMST0_true, RMST1_true, deltaRMST_true) {
  
  na <- NA_real_
  
  data.frame(
    iter = iter, seed = seed, method = method,
    match_failed = TRUE,
    
    # Cox model 0
    beta_hat_0 = na, se_hat_0 = na, ci_low_0 = na, ci_high_0 = na,
    cover_0 = na, ic_len_0 = na, p_val_0 = na, sign_pval_0 = na, bias_beta_0 = na,
    # Cox model 1
    beta_hat_1 = na, se_hat_1 = na, ci_low_1 = na, ci_high_1 = na,
    cover_1 = na, ic_len_1 = na, p_val_1 = na, sign_pval_1 = na, bias_beta_1 = na,
    # Cox model 2
    beta_hat_2 = na, se_hat_2 = na, ci_low_2 = na, ci_high_2 = na,
    cover_2 = na, ic_len_2 = na, p_val_2 = na, sign_pval_2 = na, bias_beta_2 = na,
    
    # Endpoint estimates
    trial_z1 = na, trial_z0 = na, matched = na, registry = na,
    deltaS_obs = na, deltaRMST_obs = na,
    deltaS_registry_obs = na, deltaRMST_registry_obs = na,
    deltaS_oracle_obs = na, deltaRMST_oracle_obs = na,
    
    # Latent endpoint values
    trial_true_z1 = na, trial_true_z0 = na,
    matched_true_z0 = na, registry_true = na, n_atrisk_trial = na, 
    
    # Truth (questi NON sono NA: sono parametri noti dello scenario)
    tau = tau,
    beta_unadj_star = beta_unadj_star,
    beta_adj_prog_star = beta_adj_prog_star,
    beta_adj_prop_star = beta_adj_prop_star,
    S0_true_tau = S0_true_tau, S1_true_tau = S1_true_tau,
    deltaS_true = deltaS_true,
    RMST0_true = RMST0_true, RMST1_true = RMST1_true, deltaRMST_true = deltaRMST_true,
    
    # Survival biases
    bias_S0_matched = na, bias_S0_registry = na, bias_S1_trial = na,
    bias_deltaS_matched = na, bias_deltaS_registry = na,
    
    # RMST estimates
    rmst_trial_z1 = na, rmst_trial_z0 = na, rmst_matched = na, rmst_registry = na,
    rmst_trial_true_z1 = na, rmst_trial_true_z0 = na,
    rmst_matched_true_z0 = na, rmst_registry_true = na,
    
    # RMST biases
    bias_RMST0_matched = na, bias_RMST0_registry = na, bias_RMST1_trial = na,
    bias_deltaRMST_matched = na, bias_deltaRMST_registry = na,
    
    # Matching error
    err_S0_matched_vs_trialcf = na, err_S0_trial_vs_registry = na,
    err_RMST0_matched_vs_trialcf = na, err_RMST0_trial_vs_registry = na,
    
    # Bias decomposition survival
    err_trial_z1_est = na, err_matched_z0_est = na, err_registry_z0_est = na,
    err_matching_z0 = na, err_registry_selection_z0 = na,
    sampling_error_deltaS = na,
    bias_deltaS_decomp_matched = na, bias_deltaS_decomp_naive = na,
    
    # Bias decomposition RMST
    err_trial_z1_rmst_est = na, err_matched_z0_rmst_est = na,
    err_registry_z0_rmst_est = na, err_matching_rmst_z0 = na,
    err_registry_selection_rmst = na, sampling_error_deltaRMST = na,
    bias_deltaRMST_decomp_matched = na, bias_deltaRMST_decomp_naive = na,
    
    # Balance
    score_smd_before = na, score_smd_matched = na,
    score_smd_reduction_abs = na, score_smd_reduction_rel = na,
    score_var_ratio_before = na, score_var_ratio_matched = na,
    score_pair_distance = na,
    n_control_all = na, n_trial_all = na, n_control_matched = na,
    n_control_ess = na, n_control_discarded = na, prop_control_discarded = na,
    smd_cov_mean_matched = na, smd_cov_max_matched = na, smd_cov_min_matched = na,
    smd_cov_prop_gt_0.1 = na, smd_cov_prop_gt_0.2 = na
  )
}