compute_truth_latent <- function(target_pop, tau = 12) {
  
  # exact survival truth from latent event times
  S0_tau_true <- mean(target_pop$T_true_z0 > tau, na.rm = TRUE)
  S1_tau_true <- mean(target_pop$T_true_z1 > tau, na.rm = TRUE)
  
  # exact RMST truth from latent event times
  RMST0_true <- mean(pmin(target_pop$T_true_z0, tau), na.rm = TRUE)
  RMST1_true <- mean(pmin(target_pop$T_true_z1, tau), na.rm = TRUE)
  
  
  list(
    tau = tau,
    S0_tau_true = S0_tau_true,
    S1_tau_true = S1_tau_true,
    deltaS_tau_true = S1_tau_true - S0_tau_true,
    RMST0_true = RMST0_true,
    RMST1_true = RMST1_true,
    deltaRMST_true = RMST1_true - RMST0_true
  )
}



compute_truth_superpop <- function(registry,
                                   trial_target,
                                   covs,
                                   tau = 12) {
  
  # Registry under control
  reg0 <- registry %>%
    mutate(
      Z = 0,
      S = 0,
      time = time_z0,
      event = event_z0
    )
  
  # Trial target under treatment
  trial1 <- trial_target %>%
    mutate(
      Z = 1,
      S = 1,
      time = time_z1,
      event = event_z1
    )
  
  # Trial target under control counterfactual
  trial0 <- trial_target %>%
    mutate(
      Z = 0,
      S = 1,
      time = time_z0,
      event = event_z0
    )
  
  # Dataset for Cox treatment-effect truth
  data_full <- bind_rows(
    reg0,
    trial1
  )
  
  # -----------------------------------------
  # 1. True survival and RMST in trial target
  # -----------------------------------------
  S0_true_tau <- mean(trial_target$T_true_z0 > tau, na.rm = TRUE)
  S1_true_tau <- mean(trial_target$T_true_z1 > tau, na.rm = TRUE)
  
  RMST0_true <- mean(pmin(trial_target$T_true_z0, tau), na.rm = TRUE)
  RMST1_true <- mean(pmin(trial_target$T_true_z1, tau), na.rm = TRUE)
  
  deltaS_true <- S1_true_tau - S0_true_tau
  deltaRMST_true <- RMST1_true - RMST0_true
  
  median_surv0_tau <- median(trial_target$T_true_z0, na.rm=TRUE)
  median_surv1_tau <- median(trial_target$T_true_z1, na.rm=TRUE)
  # -----------------------------------------
  # 2. Marginal / unadjusted Cox truth
  # -----------------------------------------
  fit_unadj <- coxph(
    Surv(time, event) ~ Z,
    data = data_full
  )
  
  beta_unadj_star <- unname(coef(fit_unadj)["Z"])
  
  # -----------------------------------------
  # 3. Prognostic-score-adjusted Cox truth
  # -----------------------------------------
  prog_formula <- as.formula(
    paste0("Surv(time_z0, event_z0) ~ ", paste(covs, collapse = " + "))
  )
  
  fit_prog <- coxph(
    prog_formula,
    data = registry
  )
  
  data_full$prog_score <- predict(
    fit_prog,
    newdata = data_full,
    type = "lp"
  )
  
  fit_adj_prog <- coxph(
    Surv(time, event) ~ Z + prog_score,
    data = data_full
  )
  
  beta_adj_prog_star <- unname(coef(fit_adj_prog)["Z"])
  
  # -----------------------------------------
  # 4. Propensity-score-adjusted Cox truth
  # -----------------------------------------
  ps_formula <- as.formula(
    paste0("S ~ ", paste(covs, collapse = " + "))
  )
  
  fit_ps <- glm(
    ps_formula,
    data = data_full,
    family = binomial()
  )
  
  data_full$prop_score <- predict(
    fit_ps,
    newdata = data_full,
    type = "response"
  )
  
  fit_adj_prop <- coxph(
    Surv(time, event) ~ Z + prop_score,
    data = data_full
  )
  
  beta_adj_prop_star <- unname(coef(fit_adj_prop)["Z"])
  
  # -----------------------------------------
  # Output
  # -----------------------------------------
  tibble(
    beta_unadj_star = beta_unadj_star,
    beta_adj_prog_star = beta_adj_prog_star,
    beta_adj_prop_star = beta_adj_prop_star,
    S0_true_tau = S0_true_tau,
    S1_true_tau = S1_true_tau,
    deltaS_true = deltaS_true,
    RMST0_true = RMST0_true,
    RMST1_true = RMST1_true,
    deltaRMST_true = deltaRMST_true,
    median_surv0_tau = median_surv0_tau,
    median_surv1_tau = median_surv1_tau
  )
}
