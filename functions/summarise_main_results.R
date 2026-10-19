
## summarise_main_results — con Monte Carlo standard error (ADEMP / Morris et al. 2019)
##
## MCSE per misura (n = repliche NON mancanti per quella misura):
##   - bias  = mean(b)             -> MCSE = sd(b) / sqrt(n)
##   - rmse  = sqrt(mean(b^2))     -> MCSE = sd(b^2) / (2 * rmse * sqrt(n))   [delta method]
##   - cover = mean(I)            -> MCSE = sqrt(p(1-p)/n)
##   - reject(power/typeI) = mean(I) -> MCSE = sqrt(p(1-p)/n)
##
## NB: con match_failed le righe fallite hanno NA. Il denominatore del MCSE è il
## numero di repliche valide PER MISURA, non n(). Per questo ogni misura usa il
## proprio conteggio di non-NA invece di n_iter globale.
##
## metto la mean dela variabile come ULTIMA nel suo blocco: prima si calcolano MCSE e RMSE (che leggono il
## vettore delle repliche), poi si sovrascrive il nome con la media scalare.
##
## MCSE per misura (n = repliche NON mancanti per quella misura):
##   bias  = mean(b)            -> sd(b)/sqrt(n)
##   rmse  = sqrt(mean(b^2))    -> sd(b^2)/(2*rmse*sqrt(n))     [delta method]
##   prop  = mean(I)            -> sqrt(p(1-p)/n)               [cover, type I, power]

summarise_main_results <- function(res_boot,
                                   group_vars = c("method"),
                                   alpha = 0.10,
                                   return_flextable = FALSE) {
  
  mcse_bias <- function(b) {
    n <- sum(!is.na(b)); if (n < 2) return(NA_real_)
    sd(b, na.rm = TRUE) / sqrt(n)
  }
  mcse_rmse <- function(b) {
    n  <- sum(!is.na(b)); rm <- sqrt(mean(b^2, na.rm = TRUE))
    if (is.na(rm) || rm == 0 || n < 2) return(NA_real_)
    sd(b^2, na.rm = TRUE) / (2 * rm * sqrt(n))
  }
  mcse_prop <- function(ind) {
    ind <- ind[!is.na(ind)]; n <- length(ind); if (n == 0) return(NA_real_)
    p <- mean(ind); sqrt(p * (1 - p) / n)
  }
  rej <- function(p) as.integer(p < alpha)
  
  out <- res_boot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_iter      = dplyr::n(),
      n_valid     = sum(!is.na(deltaS_obs)),
      prop_failed = mean(match_failed, na.rm = TRUE),
      
      ## ─── PRIMARY: ΔS ──────────────────────────────────────────────────────
      # (MCSE e RMSE PRIMA, media omonima PER ULTIMA)
      mean_deltaS_matched_obs = mean(deltaS_obs,    na.rm = TRUE),
      true_deltaS=mean(deltaS_true,   na.rm = TRUE),
      bias_deltaS_matched_mcse  = mcse_bias(bias_deltaS_matched),
      rmse_deltaS_matched       = sqrt(mean(bias_deltaS_matched^2, na.rm = TRUE)),
      rmse_deltaS_matched_mcse  = mcse_rmse(bias_deltaS_matched),
      total_bias_deltaS_matched        = mean(bias_deltaS_matched,       na.rm = TRUE),
      trial_est_error_S                = mean(err_trial_z1_est,          na.rm = TRUE),
      matched_control_est_error_S      = mean(err_matched_z0_est,        na.rm = TRUE),
      matching_error_S                 = mean(err_matching_z0,           na.rm = TRUE),
      sampling_error_deltaS            = mean(sampling_error_deltaS,     na.rm = TRUE),
      reconstructed_bias_deltaS_matched= mean(bias_deltaS_decomp_matched,na.rm = TRUE),
      
      mean_deltaS_naive_obs   = mean(deltaS_registry_obs,     na.rm = TRUE),
      bias_deltaS_registry_mcse = mcse_bias(bias_deltaS_registry),
      rmse_deltaS_registry      = sqrt(mean(bias_deltaS_registry^2, na.rm = TRUE)),
      rmse_deltaS_registry_mcse = mcse_rmse(bias_deltaS_registry),
      total_bias_deltaS_naive          = mean(bias_deltaS_registry,          na.rm = TRUE),
      registry_est_error_S             = mean(err_registry_z0_est,           na.rm = TRUE),
      registry_selection_error_S       = mean(err_registry_selection_z0,     na.rm = TRUE),
      reconstructed_bias_deltaS_naive  = mean(bias_deltaS_decomp_naive,      na.rm = TRUE),
      
      ## ─── CO-PRIMARY: ΔRMST ────────────────────────────────────────────────
      mean_deltaRMST_matched_obs            = mean((rmst_trial_z1-rmst_matched),     na.rm = TRUE),
      true_deltaRMST   = mean(deltaRMST_oracle_obs,            na.rm = TRUE),
      bias_deltaRMST_matched_mcse  = mcse_bias(bias_deltaRMST_matched),
      rmse_deltaRMST_matched       = sqrt(mean(bias_deltaRMST_matched^2, na.rm = TRUE)),
      rmse_deltaRMST_matched_mcse  = mcse_rmse(bias_deltaRMST_matched),
      total_bias_deltaRMST_matched        = mean(bias_deltaRMST_matched,        na.rm = TRUE),
      trial_est_error_RMST                = mean(err_trial_z1_rmst_est,         na.rm = TRUE),
      matched_control_est_error_RMST      = mean(err_matched_z0_rmst_est,       na.rm = TRUE),
      matching_error_RMST                 = mean(err_matching_rmst_z0,          na.rm = TRUE),
      sampling_error_deltaRMST            = mean(sampling_error_deltaRMST,      na.rm = TRUE),
      reconstructed_bias_deltaRMST_matched= mean(bias_deltaRMST_decomp_matched, na.rm = TRUE),
      
      mean_deltaRMST_naive_obs      = mean((rmst_trial_z1-rmst_registry),     na.rm = TRUE),
      bias_deltaRMST_registry_mcse = mcse_bias(bias_deltaRMST_registry),
      rmse_deltaRMST_registry      = sqrt(mean(bias_deltaRMST_registry^2, na.rm = TRUE)),
      rmse_deltaRMST_registry_mcse = mcse_rmse(bias_deltaRMST_registry),
      total_bias_deltaRMST_naive          = mean(bias_deltaRMST_registry,       na.rm = TRUE),
      registry_est_error_RMST             = mean(err_registry_z0_rmst_est,      na.rm = TRUE),
      registry_selection_error_RMST       = mean(err_registry_selection_rmst,   na.rm = TRUE),
      reconstructed_bias_deltaRMST_naive  = mean(bias_deltaRMST_decomp_naive,   na.rm = TRUE),
      
      ## ----SURVIVAL MATCHED -----------------------------------------------
      # ── Survival: latent true values ──────────────────────────────────────
      # Positivo = benchmark sovrastima il contraffattuale del trial
      # Negativo = benchmark sottostima (caso atteso in S4)
      mean_trial_true_z1            = mean(trial_true_z1,   na.rm = TRUE),
      sd_trial_true_z1              = sd(trial_true_z1,     na.rm = TRUE),
      mean_trial_true_z0            = mean(trial_true_z0,   na.rm = TRUE),
      sd_trial_true_z0              = sd(trial_true_z0,     na.rm = TRUE),
      mean_matched_true_z0          = mean(matched_true_z0, na.rm = TRUE),
      sd_matched_true_z0            = sd(matched_true_z0,   na.rm = TRUE),
      mean_registry_true            = mean(registry_true,   na.rm = TRUE),
      sd_registry_true              = sd(registry_true,     na.rm = TRUE),
      
      # Benchmark alignment: matched_true - trial_true e registry_true - trial_true
      # convenzione: positivo = benchmark piu favorevole del contraffattuale
      mean_diff_matched_bench_true  = mean(matched_true_z0 - trial_true_z0, na.rm = TRUE),
      mean_diff_naive_bench_true    = mean(registry_true   - trial_true_z0, na.rm = TRUE),
      
      # ── RMST: latent true values ──────────────────────────────────────────
      mean_rmst_trial_true_z1           = mean(rmst_trial_true_z1,   na.rm = TRUE),
      sd_rmst_trial_true_z1             = sd(rmst_trial_true_z1,     na.rm = TRUE),
      mean_rmst_trial_true_z0           = mean(rmst_trial_true_z0,   na.rm = TRUE),
      sd_rmst_trial_true_z0             = sd(rmst_trial_true_z0,     na.rm = TRUE),
      mean_rmst_matched_true_z0         = mean(rmst_matched_true_z0, na.rm = TRUE),
      sd_rmst_matched_true_z0           = sd(rmst_matched_true_z0,   na.rm = TRUE),
      mean_rmst_registry_true           = mean(rmst_registry_true,   na.rm = TRUE),
      sd_rmst_registry_true             = sd(rmst_registry_true,     na.rm = TRUE),
      
      # Stessa convenzione survival: positivo = benchmark piu favorevole
      mean_rmst_diff_matched_bench_true = mean(rmst_matched_true_z0 - rmst_trial_true_z0, na.rm = TRUE),
      mean_rmst_diff_naive_bench_true   = mean(rmst_registry_true   - rmst_trial_true_z0, na.rm = TRUE),
      
      # ── Survival: stime KM osservate ─────────────────────────────────────
      mean_trial_z1_est = mean(trial_z1, na.rm = TRUE),
      sd_trial_z1_est   = sd(trial_z1,   na.rm = TRUE),
      mean_trial_z0_est = mean(trial_z0, na.rm = TRUE),
      sd_trial_z0_est   = sd(trial_z0,   na.rm = TRUE),
      mean_matched_est  = mean(matched,  na.rm = TRUE),
      sd_matched_est    = sd(matched,    na.rm = TRUE),
      mean_registry_est = mean(registry, na.rm = TRUE),
      sd_registry_est   = sd(registry,   na.rm = TRUE),
      
      # ── RMST: stime osservate ─────────────────────────────────────────────
      mean_rmst_trial_z1 = mean(rmst_trial_z1, na.rm = TRUE),
      sd_rmst_trial_z1   = sd(rmst_trial_z1,   na.rm = TRUE),
      mean_rmst_trial_z0 = mean(rmst_trial_z0, na.rm = TRUE),
      sd_rmst_trial_z0   = sd(rmst_trial_z0,   na.rm = TRUE),
      mean_rmst_matched  = mean(rmst_matched,  na.rm = TRUE),
      sd_rmst_matched    = sd(rmst_matched,    na.rm = TRUE),
      mean_rmst_registry = mean(rmst_registry, na.rm = TRUE),
      sd_rmst_registry   = sd(rmst_registry,   na.rm = TRUE),
      
      
      
      ## ─── OPERATING CHARACTERISTICS (nomi nuovi: nessun shadowing) ──────────
      reject_one_sample_matched       = mean(rej(one_sample_lr_p_one_sided_better_matched),  na.rm = TRUE),
      reject_one_sample_matched_mcse  = mcse_prop(rej(one_sample_lr_p_one_sided_better_matched)),
      reject_one_sample_registry      = mean(rej(one_sample_lr_p_one_sided_better_registry), na.rm = TRUE),
      reject_one_sample_registry_mcse = mcse_prop(rej(one_sample_lr_p_one_sided_better_registry)),
      reject_two_sample_matched       = mean(rej(two_sample_lr_p_one_sided_matched), na.rm = TRUE),
      reject_two_sample_matched_mcse  = mcse_prop(rej(two_sample_lr_p_one_sided_matched)),
      reject_two_sample_naive         = mean(rej(two_sample_lr_p_one_sided_naive),   na.rm = TRUE),
      reject_two_sample_naive_mcse    = mcse_prop(rej(two_sample_lr_p_one_sided_naive)),
      
      
      ## ─── Medie assolute stime benchmark (per build_table6a) ───────────────
      mean_trial_z1_est    = mean(trial_z1,            na.rm = TRUE),
      mean_trial_z0_est    = mean(trial_z0,            na.rm = TRUE),
      mean_matched_est     = mean(matched,             na.rm = TRUE),
      mean_registry_est    = mean(registry,            na.rm = TRUE),
      mean_rmst_trial_z1   = mean(rmst_trial_z1,       na.rm = TRUE),
      mean_rmst_trial_z0   = mean(rmst_trial_z0,       na.rm = TRUE),
      mean_rmst_matched    = mean(rmst_matched,        na.rm = TRUE),
      mean_rmst_registry   = mean(rmst_registry,       na.rm = TRUE),
      mean_deltaS_matched  = mean(deltaS_obs,          na.rm = TRUE),
      mean_deltaS_registry = mean(deltaS_registry_obs, na.rm = TRUE),
      mean_deltaS_oracle   = mean(deltaS_oracle_obs,   na.rm = TRUE),
      mean_deltaRMST_matched  = mean(deltaRMST_obs,          na.rm = TRUE),
      mean_deltaRMST_registry = mean(deltaRMST_registry_obs, na.rm = TRUE),
      mean_deltaRMST_oracle   = mean(deltaRMST_oracle_obs,   na.rm = TRUE),
      
      .groups = "drop"
    )
  
  if (return_flextable) flextable::flextable(out) else out
}






