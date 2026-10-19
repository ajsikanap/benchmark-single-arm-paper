one_resample2 <- function(iter,
                          seeds,
                          trial_pool,
                          B = "all",   # one among: "all", "1", "0", "unfiltered"
                          registry,
                          prog_model,
                          n_trial,
                          covs,
                          truth,
                          ratio_controls = 1,
                          tau,
                          effect_setting = NULL) {  # ignorato: ora calcola SEMPRE entrambi
  
  set.seed(seeds[iter])
  
  # Ensure truth is a single scenario row
  if (nrow(truth) != 1) {
    stop("truth must contain exactly one scenario-specific row.")
  }
  
  # Extract truth values as scalars
  beta_unadj_star    <- truth$beta_unadj_star
  beta_adj_prog_star <- truth$beta_adj_prog_star
  beta_adj_prop_star <- truth$beta_adj_prop_star
  
  S0_true_tau     <- truth$S0_true_tau
  S1_true_tau     <- truth$S1_true_tau
  deltaS_true     <- truth$deltaS_true
  
  RMST0_true      <- truth$RMST0_true
  RMST1_true      <- truth$RMST1_true
  deltaRMST_true  <- truth$deltaRMST_true
  
  # ── Campiona il trial UNA volta sola ────────────────────────────────────────
  # trial_b contiene SIA z0 (null) SIA z1 (alternative): time_z0/event_z0 e
  # time_z1/event_z1. Lo stesso campione serve per entrambi i setting, così
  # null e alt condividono trial e benchmark (varianza Monte Carlo ridotta).
  trial_b <- dplyr::slice_sample(trial_pool, n = n_trial)
  
  # ── Helper: valuta i due log-rank dato l'outcome del trial vs un benchmark ──
  # Restituisce i p-value one-sided (treatment beneficial -> meno eventi).
  # <<< MODIFICATO: ora restituisce anche O, E, V, Z del one-sample log-rank.
  #     one_sample_logrank_km() li calcolava gia'; venivano scartati qui.
  #     Servono per il tipping point:
  #         theta* = (z_alpha + sqrt(z_alpha^2 + 4*O))^2 / (4*E)
  #         delta* = log(theta*)
  #     Senza O ed E, delta* NON e' ricostruibile a posteriori dal so
  #     
  eval_tests <- function(tr_time, tr_event,
                         hist_time, hist_event,    # benchmark per one-sample
                         bench_df) {                # data.frame con time/event/source per two-sample
    
    # one-sample log-rank (benchmark trattato come fisso)
    one_lr <- one_sample_logrank_km(
      time = tr_time, event = tr_event,
      hist_time = hist_time, hist_event = hist_event
    )
    
    # two-sample log-rank one-sided
    two_df <- rbind(
      data.frame(time = tr_time, event = tr_event, source = 1),
      bench_df |> dplyr::transmute(time, event, source = 0)
    )
    lr2 <- survival::survdiff(survival::Surv(time, event) ~ source, data = two_df)
    O_2s <- lr2$obs; E_2s <- lr2$exp; V_2s <- lr2$var
    z_2s <- (O_2s[2] - E_2s[2]) / sqrt(V_2s[2, 2])
    two_p <- pnorm(z_2s, lower.tail = TRUE)
    
    list(
      one_sample_p = one_lr$one_sample_lr_p_one_sided_better,
      two_sample_p = two_p,
      one_sample_O = one_lr$one_sample_lr_O,
      one_sample_E = one_lr$one_sample_lr_E,
      one_sample_V = one_lr$one_sample_lr_V,
      one_sample_Z = one_lr$one_sample_lr_Z
    )
  }
  
  # aggiunto: diagnostica del trial, PER SETTING.
  #     O dipende dall'outcome (event_z0 vs event_z1) ed E dai tempi di
  #     follow-up (time_z0 vs time_z1): entrambi differiscono tra null e alt.
  #     Quindi NON vanno in `common`.
  n_events_trial_null   <- sum(trial_b$event_z0)
  n_events_trial_alt    <- sum(trial_b$event_z1)
  fup_median_trial_null <- stats::median(trial_b$time_z0)
  fup_median_trial_alt  <- stats::median(trial_b$time_z1)
  n_atrisk_trial        <- nrow(trial_b)
  
  # ── Loop sui metodi: matching UNA volta, test DUE volte (null + alt) ─────────
  res_est <- lapply(c("prog", "prop"), function(method) {
    
    # ---- 1. Scoring (indipendente dall'outcome osservato) --------------------
    match_data <- score_data(
      trial    = trial_b,
      registry = registry,
      prog_model = prog_model
    )
    
    # ---- 2. Matching (indipendente dall'outcome osservato) -------------------
    matched_list <- tryCatch(
      run_matching(
        B = B,
        match_data = match_data,
        covariates = covs,
        method = method,
        ratio_controls = ratio_controls
      ),
      error = function(e) {
        warning(sprintf("matchit failed (iter=%d, method=%s): %s",
                        iter, method, conditionMessage(e)))
        NULL
      }
    )
    
    match_failed <- is.null(matched_list) ||
      is.null(matched_list$matched) ||
      nrow(matched_list$matched) == 0
    
    # ---- 2b. Se il matching fallisce: restituisci due righe NA (null + alt) ---
    if (match_failed) {
      na_row <- make_na_row(
        iter = iter, seed = seeds[iter], method = method, tau = tau,
        beta_unadj_star = beta_unadj_star,
        beta_adj_prog_star = beta_adj_prog_star,
        beta_adj_prop_star = beta_adj_prop_star,
        S0_true_tau = S0_true_tau, S1_true_tau = S1_true_tau,
        deltaS_true = deltaS_true,
        RMST0_true = RMST0_true, RMST1_true = RMST1_true,
        deltaRMST_true = deltaRMST_true
      )
      na_row$one_sample_lr_p_one_sided_better_matched  <- NA_real_
      na_row$two_sample_lr_p_one_sided_matched         <- NA_real_
      na_row$one_sample_lr_p_one_sided_better_registry <- NA_real_
      na_row$two_sample_lr_p_one_sided_naive           <- NA_real_
    
      
      # Aggiunto:: le stesse colonne devono esistere anche nella riga NA, altrimenti
      #     il bind_rows finale produce ncol diversi
      #     
      na_row$one_sample_lr_O_matched  <- NA_real_
      na_row$one_sample_lr_E_matched  <- NA_real_
      na_row$one_sample_lr_V_matched  <- NA_real_
      na_row$one_sample_lr_Z_matched  <- NA_real_
      na_row$one_sample_lr_O_registry <- NA_real_
      na_row$one_sample_lr_E_registry <- NA_real_
      na_row$one_sample_lr_V_registry <- NA_real_
      na_row$one_sample_lr_Z_registry <- NA_real_
      na_row$n_events_trial           <- NA_real_
      na_row$n_atrisk_trial           <- NA_real_
      na_row$fup_median_trial         <- NA_real_
      
      
      
      return(dplyr::bind_rows(
        dplyr::mutate(na_row, effect_setting = "null"),
        dplyr::mutate(na_row, effect_setting = "alternative")
      ))
    }
    
    matched_controls <- matched_list$matched |> dplyr::filter(source == 0)
    
    # ---- 3. Endpoint assessment (calcolato sotto z1 = "alternative") ---------
    # assess_endpoint usa trial$time/event per il braccio trattato: impostiamo
    # l'outcome alternative (z1). I valori latenti (trial_true_z0/z1, ecc.) e i
    # benchmark NON dipendono dal setting, quindi la decomposizione del bias
    # vale per "alternative" come richiesto.
    trial_alt <- trial_b |>
      dplyr::mutate(time = time_z1, event = event_z1)
    
    surv <- assess_endpoint(
      trial    = trial_alt,
      matched  = matched_list$matched,
      registry = registry,
      tau      = tau
    )
    
    # ---- 4. Cox estimates (sotto alternative) --------------------------------
    est1 <- estimate_treatment(matched_list$matched, truth$beta_unadj_star,    "source_only")
    est2 <- estimate_treatment(matched_list$matched, truth$beta_adj_prog_star, "source_prog_score")
    est3 <- estimate_treatment(
      matched_list$matched,
      beta_true = if (method == "prop") truth$beta_adj_prop_star else truth$beta_adj_prog_star,
      method = "source_score"
    )
    
    # ---- 5. TEST: valutati per ENTRAMBI i setting sullo stesso benchmark ------
    # NULL: trial usa z0  |  ALTERNATIVE: trial usa z1
    tests_null <- eval_tests(
      tr_time = trial_b$time_z0, tr_event = trial_b$event_z0,
      hist_time = matched_controls$time, hist_event = matched_controls$event,
      bench_df = matched_controls
    )
    tests_alt <- eval_tests(
      tr_time = trial_b$time_z1, tr_event = trial_b$event_z1,
      hist_time = matched_controls$time, hist_event = matched_controls$event,
      bench_df = matched_controls
    )
    
    # one-sample / two-sample contro il REGISTRO pieno (naive), per setting
    tests_null_reg <- eval_tests(
      tr_time = trial_b$time_z0, tr_event = trial_b$event_z0,
      hist_time = registry$time, hist_event = registry$event,
      bench_df = registry
    )
    tests_alt_reg <- eval_tests(
      tr_time = trial_b$time_z1, tr_event = trial_b$event_z1,
      hist_time = registry$time, hist_event = registry$event,
      bench_df = registry
    )
    
    # ---- 6. Osservati / bias (alternative) -----------------------------------
    deltaS_obs            <- surv$trial_z1 - surv$matched_z0
    deltaRMST_obs         <- surv$rmst_trial_z1 - surv$rmst_matched_z0
    deltaS_registry_obs   <- surv$trial_z1 - surv$registry
    deltaRMST_registry_obs<- surv$rmst_trial_z1 - surv$rmst_registry
    deltaS_oracle_obs     <- surv$trial_z1 - surv$trial_true_z0
    deltaRMST_oracle_obs  <- surv$rmst_trial_z1 - surv$rmst_trial_true_z0
    
    cov_names <- intersect(names(covs), names(match_data))
    if (length(cov_names) > 0) {
      smd_dt <- cobalt::bal.tab(
        matched_list$m.out,
        data = matched_list$filtered_matched,
        addl = cov_names, un = FALSE
      )
    } else {
      smd_dt <- NULL
    }
    
    # Cox biases
    bias_beta_0 <- est1$beta_hat - beta_unadj_star
    bias_beta_1 <- est2$beta_hat - beta_adj_prog_star
    bias_beta_2 <- est3$beta_hat -
      if (method == "prop") beta_adj_prop_star else beta_adj_prog_star
    
    # Matching error on survival scale
    err_S0_matched_vs_trialcf <- surv$matched_true_z0 - surv$trial_true_z0
    err_S0_trial_vs_registry  <- surv$trial_true_z0   - surv$registry_true
    
    # Survival bias vs target truth
    bias_S0_matched  <- surv$matched_z0 - truth$S0_true_tau
    bias_S0_registry <- surv$registry   - truth$S0_true_tau
    bias_S1_trial    <- surv$trial_z1   - truth$S1_true_tau
    
    bias_deltaS_matched  <- (surv$trial_z1 - surv$matched_z0) - truth$deltaS_true
    bias_deltaS_registry <- (surv$trial_z1 - surv$registry)   - truth$deltaS_true
    
    sampling_error_deltaS <- (surv$trial_true_z1 - surv$trial_true_z0) - truth$deltaS_true
    
    # Survival decomposition
    err_trial_z1_est          <- surv$trial_z1   - surv$trial_true_z1
    err_matched_z0_est        <- surv$matched_z0 - surv$matched_true_z0
    err_registry_z0_est       <- surv$registry   - surv$registry_true
    err_matching_z0           <- surv$matched_true_z0 - surv$trial_true_z0
    err_registry_selection_z0 <- surv$registry_true   - surv$trial_true_z0 # questo e' sbagliato, dovrebbe essere il contratio
    
    bias_deltaS_decomp_matched <- err_trial_z1_est - err_matched_z0_est -
      err_matching_z0 + sampling_error_deltaS
    bias_deltaS_decomp_naive   <- err_trial_z1_est - err_registry_z0_est -
      err_registry_selection_z0 + sampling_error_deltaS
    
    # RMST availability
    has_rmst <- all(c("rmst_trial_true_z1", "rmst_trial_true_z0",
                      "rmst_matched_true_z0", "rmst_registry_true") %in% names(surv))
    
    if (has_rmst) {
      bias_RMST0_matched  <- surv$rmst_matched_z0 - truth$RMST0_true
      bias_RMST0_registry <- surv$rmst_registry   - truth$RMST0_true
      bias_RMST1_trial    <- surv$rmst_trial_z1   - truth$RMST1_true
      
      bias_deltaRMST_matched  <- (surv$rmst_trial_z1 - surv$rmst_matched_z0) - truth$deltaRMST_true
      bias_deltaRMST_registry <- (surv$rmst_trial_z1 - surv$rmst_registry)   - truth$deltaRMST_true
      
      err_trial_z1_rmst_est       <- surv$rmst_trial_z1      - surv$rmst_trial_true_z1
      err_matched_z0_rmst_est     <- surv$rmst_matched_z0    - surv$rmst_matched_true_z0
      err_registry_z0_rmst_est    <- surv$rmst_registry      - surv$rmst_registry_true
      err_matching_rmst_z0        <- surv$rmst_matched_true_z0 - surv$rmst_trial_true_z0
      err_registry_selection_rmst <- surv$rmst_registry_true   - surv$rmst_trial_true_z0
      sampling_error_deltaRMST    <- (surv$rmst_trial_true_z1 - surv$rmst_trial_true_z0) - truth$deltaRMST_true
      
      err_RMST0_matched_vs_trialcf <- surv$rmst_matched_true_z0 - surv$rmst_trial_true_z0
      err_RMST0_trial_vs_registry  <- surv$rmst_trial_true_z0   - surv$rmst_registry_true
      
      bias_deltaRMST_decomp_matched <- err_trial_z1_rmst_est - err_matched_z0_rmst_est -
        err_matching_rmst_z0 + sampling_error_deltaRMST
      bias_deltaRMST_decomp_naive   <- err_trial_z1_rmst_est - err_registry_z0_rmst_est -
        err_registry_selection_rmst + sampling_error_deltaRMST
    } else {
      bias_RMST0_matched <- bias_RMST0_registry <- bias_RMST1_trial <-
        bias_deltaRMST_matched <- bias_deltaRMST_registry <-
        err_trial_z1_rmst_est <- err_matched_z0_rmst_est <- err_registry_z0_rmst_est <-
        err_matching_rmst_z0 <- err_registry_selection_rmst <- sampling_error_deltaRMST <-
        err_RMST0_matched_vs_trialcf <- err_RMST0_trial_vs_registry <-
        bias_deltaRMST_decomp_matched <- bias_deltaRMST_decomp_naive <- NA_real_
    }
    
    match_diag <- extract_match_diagnostics(
      m.out = matched_list$m.out, smd_dt = smd_dt
    )
    
    # ── Blocco "comune": tutto ciò che NON dipende dal setting ────────────────
    # (bias, decomposizione, diagnostica matching, verità).
    common <- data.frame(
      iter   = iter,
      seed   = seeds[iter],
      method = method,
      match_failed = FALSE,
      
      # Cox model 0
      beta_hat_0 = est1$beta_hat, se_hat_0 = est1$se_hat,
      ci_low_0 = est1$ci_low, ci_high_0 = est1$ci_high,
      cover_0 = est1$cover, ic_len_0 = est1$ic_len,
      p_val_0 = est1$p_val, sign_pval_0 = est1$sign_pval, bias_beta_0 = bias_beta_0,
      # Cox model 1
      beta_hat_1 = est2$beta_hat, se_hat_1 = est2$se_hat,
      ci_low_1 = est2$ci_low, ci_high_1 = est2$ci_high,
      cover_1 = est2$cover, ic_len_1 = est2$ic_len,
      p_val_1 = est2$p_val, sign_pval_1 = est2$sign_pval, bias_beta_1 = bias_beta_1,
      # Cox model 2
      beta_hat_2 = est3$beta_hat, se_hat_2 = est3$se_hat,
      ci_low_2 = est3$ci_low, ci_high_2 = est3$ci_high,
      cover_2 = est3$cover, ic_len_2 = est3$ic_len,
      p_val_2 = est3$p_val, sign_pval_2 = est3$sign_pval, bias_beta_2 = bias_beta_2,
      
      # Endpoint estimates (alternative)
      trial_z1 = surv$trial_z1,
      trial_z0 = surv$trial_z0,
      matched = surv$matched_z0,
      registry = surv$registry,
      deltaS_obs = deltaS_obs, 
      deltaRMST_obs = deltaRMST_obs,
      deltaS_registry_obs = deltaS_registry_obs,
      deltaRMST_registry_obs = deltaRMST_registry_obs,
      deltaS_oracle_obs = deltaS_oracle_obs,
      deltaRMST_oracle_obs = deltaRMST_oracle_obs,
      
      # Latent endpoint values
      trial_true_z1 = surv$trial_true_z1, trial_true_z0 = surv$trial_true_z0,
      matched_true_z0 = surv$matched_true_z0, registry_true = surv$registry_true,
      
      # <<< NUOVO: numerosita' del trial (non dipende dal setting)
      n_atrisk_trial = n_atrisk_trial,
      
      # Truth
      tau = tau,
      beta_unadj_star = beta_unadj_star,
      beta_adj_prog_star = beta_adj_prog_star,
      beta_adj_prop_star = beta_adj_prop_star,
      S0_true_tau = S0_true_tau, S1_true_tau = S1_true_tau,
      deltaS_true = deltaS_true,
      RMST0_true = RMST0_true, RMST1_true = RMST1_true, 
      deltaRMST_true = deltaRMST_true,
      
      # Survival biases
      bias_S0_matched = bias_S0_matched, 
      bias_S0_registry = bias_S0_registry,
      bias_S1_trial = bias_S1_trial,
      bias_deltaS_matched = bias_deltaS_matched,
      bias_deltaS_registry = bias_deltaS_registry,
      
      # RMST estimates
      rmst_trial_z1 = surv$rmst_trial_z1, rmst_trial_z0 = surv$rmst_trial_z0,
      rmst_matched = surv$rmst_matched_z0, rmst_registry = surv$rmst_registry,
      rmst_trial_true_z1 = surv$rmst_trial_true_z1, rmst_trial_true_z0 = surv$rmst_trial_true_z0,
      rmst_matched_true_z0 = surv$rmst_matched_true_z0, rmst_registry_true = surv$rmst_registry_true,
      
      # RMST biases
      bias_RMST0_matched = bias_RMST0_matched, bias_RMST0_registry = bias_RMST0_registry,
      bias_RMST1_trial = bias_RMST1_trial,
      bias_deltaRMST_matched = bias_deltaRMST_matched, bias_deltaRMST_registry = bias_deltaRMST_registry,
      
      # Matching error
      err_S0_matched_vs_trialcf = err_S0_matched_vs_trialcf,
      err_S0_trial_vs_registry = err_S0_trial_vs_registry,
      err_RMST0_matched_vs_trialcf = err_RMST0_matched_vs_trialcf,
      err_RMST0_trial_vs_registry = err_RMST0_trial_vs_registry,
      
      # Bias decomposition survival
      err_trial_z1_est = err_trial_z1_est,
      err_matched_z0_est = err_matched_z0_est,
      err_registry_z0_est = err_registry_z0_est,
      err_matching_z0 = err_matching_z0,
      err_registry_selection_z0 = err_registry_selection_z0,
      sampling_error_deltaS = sampling_error_deltaS,
      bias_deltaS_decomp_matched = bias_deltaS_decomp_matched,
      bias_deltaS_decomp_naive = bias_deltaS_decomp_naive,
      
      # Bias decomposition RMST
      err_trial_z1_rmst_est = err_trial_z1_rmst_est,
      err_matched_z0_rmst_est = err_matched_z0_rmst_est,
      err_registry_z0_rmst_est = err_registry_z0_rmst_est,
      err_matching_rmst_z0 = err_matching_rmst_z0,
      err_registry_selection_rmst = err_registry_selection_rmst,
      sampling_error_deltaRMST = sampling_error_deltaRMST,
      bias_deltaRMST_decomp_matched = bias_deltaRMST_decomp_matched,
      bias_deltaRMST_decomp_naive = bias_deltaRMST_decomp_naive,
      
      # Balance
      score_smd_before = match_diag$score_smd_before,
      score_smd_matched = match_diag$score_smd_matched,
      score_smd_reduction_abs = match_diag$score_smd_reduction_abs,
      score_smd_reduction_rel = match_diag$score_smd_reduction_rel,
      score_var_ratio_before = match_diag$score_var_ratio_before,
      score_var_ratio_matched = match_diag$score_var_ratio_matched,
      score_pair_distance = match_diag$score_pair_distance,
      n_control_all = match_diag$n_control_all,
      n_trial_all = match_diag$n_trial_all,
      n_control_matched = match_diag$n_control_matched,
      n_control_ess = match_diag$n_control_ess,
      n_control_discarded = match_diag$n_control_discarded,
      prop_control_discarded = match_diag$prop_control_discarded,
      smd_cov_mean_matched = match_diag$smd_cov_mean_matched,
      smd_cov_max_matched = match_diag$smd_cov_max_matched,
      smd_cov_min_matched = match_diag$smd_cov_min_matched,
      smd_cov_prop_gt_0.1 = match_diag$smd_cov_prop_gt_0.1,
      smd_cov_prop_gt_0.2 = match_diag$smd_cov_prop_gt_0.2
    )
    
    # ── Riga NULL: test sotto z0; bias = NA (richiesta: bias solo alternative) ─
    row_null <- common
    row_null$effect_setting                            <- "null"
    row_null$one_sample_lr_p_one_sided_better_matched  <- tests_null$one_sample_p
    row_null$two_sample_lr_p_one_sided_matched         <- tests_null$two_sample_p
    row_null$one_sample_lr_p_one_sided_better_registry <- tests_null_reg$one_sample_p
    row_null$two_sample_lr_p_one_sided_naive           <- tests_null_reg$two_sample_p
    
    
    # <<< NUOVO: componenti one-sample sotto z0
    row_null$one_sample_lr_O_matched  <- tests_null$one_sample_O
    row_null$one_sample_lr_E_matched  <- tests_null$one_sample_E
    row_null$one_sample_lr_V_matched  <- tests_null$one_sample_V
    row_null$one_sample_lr_Z_matched  <- tests_null$one_sample_Z
    row_null$one_sample_lr_O_registry <- tests_null_reg$one_sample_O
    row_null$one_sample_lr_E_registry <- tests_null_reg$one_sample_E
    row_null$one_sample_lr_V_registry <- tests_null_reg$one_sample_V
    row_null$one_sample_lr_Z_registry <- tests_null_reg$one_sample_Z
    row_null$n_events_trial           <- n_events_trial_null
    row_null$fup_median_trial         <- fup_median_trial_null
    
    # ── Riga ALTERNATIVE: test sotto z1; bias popolati ─────────────────────────
    row_alt <- common
    row_alt$effect_setting                            <- "alternative"
    row_alt$one_sample_lr_p_one_sided_better_matched  <- tests_alt$one_sample_p
    row_alt$two_sample_lr_p_one_sided_matched         <- tests_alt$two_sample_p
    row_alt$one_sample_lr_p_one_sided_better_registry <- tests_alt_reg$one_sample_p
    row_alt$two_sample_lr_p_one_sided_naive           <- tests_alt_reg$two_sample_p
    
    
    #<<< NUOVO: componenti one-sample sotto z1
    row_alt$one_sample_lr_O_matched  <- tests_alt$one_sample_O
    row_alt$one_sample_lr_E_matched  <- tests_alt$one_sample_E
    row_alt$one_sample_lr_V_matched  <- tests_alt$one_sample_V
    row_alt$one_sample_lr_Z_matched  <- tests_alt$one_sample_Z
    row_alt$one_sample_lr_O_registry <- tests_alt_reg$one_sample_O
    row_alt$one_sample_lr_E_registry <- tests_alt_reg$one_sample_E
    row_alt$one_sample_lr_V_registry <- tests_alt_reg$one_sample_V
    row_alt$one_sample_lr_Z_registry <- tests_alt_reg$one_sample_Z
    row_alt$n_events_trial           <- n_events_trial_alt
    row_alt$fup_median_trial         <- fup_median_trial_alt
    
    dplyr::bind_rows(row_null, row_alt)
  })
  
  # only temporary
 # dplyr::bind_rows(res_est)
  ncols <- vapply(res_est, ncol, integer(1))
  if (length(unique(ncols)) > 1) {
    saveRDS(list(iter = iter, seed = seeds[iter], res_est = res_est),
            sprintf("debug_iter_%d_%s.rds", iter, B))
    return(dplyr::bind_rows(res_est))   # union-bind: il run prosegue
  }
  dplyr::bind_rows(res_est)
}
