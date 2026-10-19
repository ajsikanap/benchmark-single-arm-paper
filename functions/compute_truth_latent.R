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

