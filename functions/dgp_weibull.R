# =============================================================================
# R/dgp_weibull.R
#
# Weibull proportional-hazards data-generating process (DGP) for overall
# survival in recurrent glioblastoma. Shared by the main analysis
# (analysis/main/01_dgp_weibull.Rmd) and the sensitivity analysis
# (analysis/sensitivity/S5_01_dgp.Rmd): any change here affects both.
#
# Model (k = Weibull shape, lambda0 = Weibull scale):
#   S(t | Z, X) = exp{ -lambda0 * t^k * exp(beta_Z * Z + beta' X) }
#   T = [ -log(U) / (lambda0 * exp(beta_Z * Z + beta' X)) ]^(1/k),  U ~ U(0,1)
#
# lambda0 and k are calibrated so that the marginal control survival matches
# a target median and a target S0(tau); beta_Z is calibrated so that the
# marginal treated survival matches a target median.
#
# DGP parameters are read from one row of config/simulation_structure.xlsx,
# checked with check_dgp_params(), and passed as explicit arguments.
# Calibration targets, population sizes and censoring bounds have no default
# values: they must be supplied from config/parameters.R.
# =============================================================================


# -----------------------------------------------------------------------------
# Parameter check
# -----------------------------------------------------------------------------

# Check that all DGP parameters are present, uniquely named and numeric
#
# @param pars One-row data frame of parameters read from the workbook.
# @return `pars`, invisibly. Stops with an informative error otherwise.
check_dgp_params <- function(pars) {
  required <- c(
    "p_M", "p_B", "p_sexF", "p_WHOm2", "p_multifocal", "p_pTERT", "p_EGFR",
    "p_g7l10", "p_unnec", "p_unmeas",
    paste0("p_resection_", 0:3), paste0("p_trt_baseline_", 0:2),
    paste0("p_loc_superficial_", 0:5),
    "mean_age", "sd_age", "mean_tumor_size", "sd_tumor_size",
    "Unmeas_mean", "Unmeas_sd", "rho_B_TS", "rho_U_TS",
    "int_OS", "beta_M_OS", "beta_B_OS", "beta_age_OS", "beta_WHOm2_OS",
    "beta_sex_OS", "beta_resection_OS", "beta_trt_baseline_OS",
    "beta_localization_OS", "beta_multifocal_OS", "beta_pTERT_OS",
    "beta_tumor_size_OS", "beta_EGFR_OS", "beta_g7l10_OS",
    "beta_Unmeas_OS", "beta_unmeas_binary"
  )
  
  if (nrow(pars) != 1) stop("Expected exactly one row of parameters.")
  
  dup <- unique(names(pars)[duplicated(names(pars))])
  if (length(dup)) stop("Duplicated parameter names: ", paste(dup, collapse = ", "))
  
  missing <- setdiff(required, names(pars))
  if (length(missing)) stop("Missing DGP parameters: ", paste(missing, collapse = ", "))
  
  bad <- required[!vapply(pars[required], is.numeric, logical(1))]
  if (length(bad)) stop("Non-numeric DGP parameters: ", paste(bad, collapse = ", "))
  
  invisible(pars)
}


# -----------------------------------------------------------------------------
# Baseline covariates
# -----------------------------------------------------------------------------

# Generate baseline covariates
#
# Covariates are generated separately from event times because the
# calibration operates on the distribution of the linear predictor.
#
# Tumour size and brain penetrance B are obtained from a bivariate normal
# latent pair with correlation rho_B_TS. B is obtained by thresholding its
# latent component, so that P(B = 1) = p_B; with rho_B_TS > 0, larger tumours
# are more likely to have B = 1.
#
# The order of the random draws is fixed and must not be changed: it
# determines the reproducibility of the superpopulations used in the paper.
#
# N Number of individuals.
# p_resection_0_1_2_3 Probabilities of the four extent-of-resection levels.
# p_trt_baseline_0_1_2 Probabilities of the three treatment-intensity levels.
# rho_B_TS Latent correlation between B and tumour size.
# rho_U_TS Latent correlation of the unmeasured binary confounder
#  (no effect: see note in the body).
# data.frame with one row per individual.
generate_covariates <- function(
    N,
    p_M,
    p_B,
    p_sexF,
    p_resection_0_1_2_3,
    p_trt_baseline_0_1_2,
    p_WHOm2,
    p_loc_superficial_0,
    p_loc_superficial_1,
    p_loc_superficial_2,
    p_loc_superficial_3,
    p_loc_superficial_4,
    p_loc_superficial_5,
    p_multifocal,
    mean_tumor_size,
    sd_tumor_size,
    p_pTERT,
    p_EGFR,
    p_g7l10,
    mean_age,
    sd_age,
    Unmeas_mean,
    Unmeas_sd,
    p_unnec,
    p_unmeas,
    rho_B_TS,
    rho_U_TS
) {
  N <- as.integer(N)
  
  rho_B_TS <- max(min(as.numeric(rho_B_TS)[1], 0.999), -0.999)
  rho_U_TS <- max(min(as.numeric(rho_U_TS)[1], 0.999), -0.999)
  
  # Latent pair (tumour size, B)
  Sigma_B_TS <- matrix(c(1, rho_B_TS, rho_B_TS, 1), nrow = 2, byrow = TRUE)
  latent     <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = Sigma_B_TS)
  
  TS_latent  <- latent[, 1]
  B_latent   <- latent[, 2]
  
  tumor_size <- pmax(mean_tumor_size + sd_tumor_size * TS_latent, 1)
  B          <- as.numeric(B_latent > qnorm(1-p_B)) 
  
  # Latent pair for the unmeasured binary confounder.
  # NOTE: this pair is drawn independently of TS_latent, so rho_U_TS has no
  # effect on the joint distribution of the confounder and tumour size.
  # rho_U_TS = 0 in all analyses of the paper. The draw is kept to preserve
  # the random number stream.
  Sigma_U_TS    <- matrix(c(1, rho_U_TS, rho_U_TS, 1), nrow = 2)
  latent_U      <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_U_TS)
  unmeas_binary <- as.numeric(latent_U[, 2] < qnorm(p_unmeas))
  
  # Arguments are evaluated left to right: do not reorder.
  data.frame(
    id             = seq_len(N),
    M              = rbinom(N, 1, p_M),
    B              = B,
    age            = rnorm(N, mean_age, sd_age),
    sex            = rbinom(N, 1, p_sexF),
    WHOm2          = rbinom(N, 1, p_WHOm2),
    resection      = sample(0:3, N, TRUE, p_resection_0_1_2_3),
    trt_type       = sample(0:2, N, TRUE, p_trt_baseline_0_1_2),
    localization_0 = rbinom(N, 1, p_loc_superficial_0),
    localization_1 = rbinom(N, 1, p_loc_superficial_1),
    localization_2 = rbinom(N, 1, p_loc_superficial_2),
    localization_3 = rbinom(N, 1, p_loc_superficial_3),
    localization_4 = rbinom(N, 1, p_loc_superficial_4),
    localization_5 = rbinom(N, 1, p_loc_superficial_5),
    multifocal     = rbinom(N, 1, p_multifocal),
    tumor_size     = tumor_size,
    pTERT          = rbinom(N, 1, p_pTERT),
    EGFR           = rbinom(N, 1, p_EGFR),
    g7l10          = rbinom(N, 1, p_g7l10),
    unmeas_binary  = unmeas_binary,
    unmeasured     = rnorm(N, Unmeas_mean, Unmeas_sd),
    unneccessary   = rbinom(N, 1, p_unnec)
  )
}


# -----------------------------------------------------------------------------
# Linear predictor
# -----------------------------------------------------------------------------

# Prognostic linear predictor for OS under control
#
# Positive coefficients increase the hazard, negative coefficients decrease
# it. Age enters centred at `age_ref`. All localization indicators share the
# same coefficient.
#
# @param dat Covariates, as returned by generate_covariates().
# @param age_ref Centring value for age (years).
# @return Numeric vector of linear predictors.
compute_lp_OS <- function(
    dat,
    int_OS,
    beta_M_OS,
    beta_B_OS,
    beta_age_OS,
    beta_WHOm2_OS,
    beta_sex_OS,
    beta_resection_OS,
    beta_trt_baseline_OS,
    beta_localization_OS,
    beta_multifocal_OS,
    beta_pTERT_OS,
    beta_tumor_size_OS,
    beta_EGFR_OS,
    beta_g7l10_OS,
    beta_Unmeas_OS,
    beta_unmeas_binary,
    age_ref = 60
) {
  int_OS +
    beta_M_OS            * dat$M +
    beta_B_OS            * dat$B +
    beta_age_OS          * (dat$age - age_ref) +
    beta_Unmeas_OS       * dat$unmeasured +
    beta_WHOm2_OS        * dat$WHOm2 +
    beta_sex_OS          * dat$sex +
    beta_resection_OS    * dat$resection +
    beta_trt_baseline_OS * dat$trt_type +
    beta_localization_OS * dat$localization_0 +
    beta_localization_OS * dat$localization_1 +
    beta_localization_OS * dat$localization_2 +
    beta_localization_OS * dat$localization_3 +
    beta_localization_OS * dat$localization_4 +
    beta_localization_OS * dat$localization_5 +
    beta_multifocal_OS   * dat$multifocal +
    beta_pTERT_OS        * dat$pTERT +
    beta_tumor_size_OS   * dat$tumor_size +
    beta_EGFR_OS         * dat$EGFR +
    beta_g7l10_OS        * dat$g7l10 +
    beta_unmeas_binary   * dat$unmeas_binary
}


# -----------------------------------------------------------------------------
# Marginal Weibull survival quantities
# -----------------------------------------------------------------------------

# Marginal survival at time t: mean of the individual survival probabilities
weibull_survival_marginal <- function(t, lp, lambda0, shape) {
  mean(exp(-lambda0 * t^shape * exp(lp)))
}

# Marginal median survival time: root of S(t) - 0.5
weibull_median_marginal <- function(lp, lambda0, shape) {
  uniroot(
    function(t) weibull_survival_marginal(t, lp, lambda0, shape) - 0.5,
    interval = c(1e-6, 500)
  )$root
}


# -----------------------------------------------------------------------------
# Calibration
# -----------------------------------------------------------------------------

# Calibrate lambda0 and k of the control arm
#
# Minimises the squared distance of the marginal control survival from 0.5 at
# the target median and from target_S0_tau at tau. The optimisation is on the
# log scale, to keep both parameters positive, and starts from the
# exponential model (k = 1, lambda0 = log(2) / median).
#
# @param lp0 Linear predictor under control.
# @return List with lambda0, shape, value of the objective and checks.
calibrate_weibull_control <- function(lp0, target_median_z0, target_S0_tau, tau) {
  objective <- function(par) {
    lambda0 <- exp(par[1])
    shape   <- exp(par[2])
    (weibull_survival_marginal(target_median_z0, lp0, lambda0, shape) - 0.5)^2 +
      (weibull_survival_marginal(tau, lp0, lambda0, shape) - target_S0_tau)^2
  }
  
  fit <- optim(
    par     = c(log(log(2) / target_median_z0), log(1)),
    fn      = objective,
    method  = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  
  lambda0 <- exp(fit$par[1])
  shape   <- exp(fit$par[2])
  
  list(
    lambda0         = lambda0,
    shape           = shape,
    objective       = fit$value,
    S0_median_check = weibull_survival_marginal(target_median_z0, lp0, lambda0, shape),
    S0_tau_check    = weibull_survival_marginal(tau, lp0, lambda0, shape),
    median_z0_check = weibull_median_marginal(lp0, lambda0, shape)
  )
}

# Calibrate the treatment log-hazard ratio
#
# Finds beta_Z such that the marginal treated survival equals 0.5 at the
# target treated median, with lambda0 and k fixed at their control values.
calibrate_beta_treatment <- function(lp0, lambda0, shape, target_median_z1,
                                     interval = c(-5, 2)) {
  objective <- function(beta_treatment_OS) {
    (weibull_survival_marginal(target_median_z1, lp0 + beta_treatment_OS,
                               lambda0, shape) - 0.5)^2
  }
  optimize(f = objective, interval = interval)$minimum
}

# Full calibration of the Weibull DGP
#
# Generates a calibration population of size N_calib, computes its linear
# predictor under control, then calibrates the control arm (lambda0, k) and
# the treatment effect (beta_Z).
#
# @return List with lambda0, weibull_shape, beta_treatment_OS and a
#   data.frame of checks (targets against achieved values).
calibrate_weibull_DGP <- function(
    N_calib,
    target_median_z0,
    target_S0_tau,
    target_median_z1,
    tau,
    p_M,
    p_B,
    p_sexF,
    p_resection_0_1_2_3,
    p_trt_baseline_0_1_2,
    p_WHOm2,
    p_loc_superficial_0,
    p_loc_superficial_1,
    p_loc_superficial_2,
    p_loc_superficial_3,
    p_loc_superficial_4,
    p_loc_superficial_5,
    p_multifocal,
    mean_tumor_size,
    sd_tumor_size,
    p_pTERT,
    p_EGFR,
    p_g7l10,
    mean_age,
    sd_age,
    Unmeas_mean,
    Unmeas_sd,
    p_unnec,
    p_unmeas,
    rho_B_TS,
    rho_U_TS,
    int_OS,
    beta_M_OS,
    beta_B_OS,
    beta_age_OS,
    beta_WHOm2_OS,
    beta_sex_OS,
    beta_resection_OS,
    beta_trt_baseline_OS,
    beta_localization_OS,
    beta_multifocal_OS,
    beta_pTERT_OS,
    beta_tumor_size_OS,
    beta_EGFR_OS,
    beta_g7l10_OS,
    beta_Unmeas_OS,
    beta_unmeas_binary
) {
  calib_pop <- generate_covariates(
    N                    = N_calib,
    p_M                  = p_M,
    p_B                  = p_B,
    p_sexF               = p_sexF,
    p_resection_0_1_2_3  = p_resection_0_1_2_3,
    p_trt_baseline_0_1_2 = p_trt_baseline_0_1_2,
    p_WHOm2              = p_WHOm2,
    p_loc_superficial_0  = p_loc_superficial_0,
    p_loc_superficial_1  = p_loc_superficial_1,
    p_loc_superficial_2  = p_loc_superficial_2,
    p_loc_superficial_3  = p_loc_superficial_3,
    p_loc_superficial_4  = p_loc_superficial_4,
    p_loc_superficial_5  = p_loc_superficial_5,
    p_multifocal         = p_multifocal,
    mean_tumor_size      = mean_tumor_size,
    sd_tumor_size        = sd_tumor_size,
    p_pTERT              = p_pTERT,
    p_EGFR               = p_EGFR,
    p_g7l10              = p_g7l10,
    mean_age             = mean_age,
    sd_age               = sd_age,
    Unmeas_mean          = Unmeas_mean,
    Unmeas_sd            = Unmeas_sd,
    p_unnec              = p_unnec,
    p_unmeas             = p_unmeas,
    rho_B_TS             = rho_B_TS,
    rho_U_TS             = rho_U_TS
  )
  
  lp0 <- compute_lp_OS(
    dat                  = calib_pop,
    int_OS               = int_OS,
    beta_M_OS            = beta_M_OS,
    beta_B_OS            = beta_B_OS,
    beta_age_OS          = beta_age_OS,
    beta_WHOm2_OS        = beta_WHOm2_OS,
    beta_sex_OS          = beta_sex_OS,
    beta_resection_OS    = beta_resection_OS,
    beta_trt_baseline_OS = beta_trt_baseline_OS,
    beta_localization_OS = beta_localization_OS,
    beta_multifocal_OS   = beta_multifocal_OS,
    beta_pTERT_OS        = beta_pTERT_OS,
    beta_tumor_size_OS   = beta_tumor_size_OS,
    beta_EGFR_OS         = beta_EGFR_OS,
    beta_g7l10_OS        = beta_g7l10_OS,
    beta_Unmeas_OS       = beta_Unmeas_OS,
    beta_unmeas_binary   = beta_unmeas_binary
  )
  
  control_cal <- calibrate_weibull_control(
    lp0              = lp0,
    target_median_z0 = target_median_z0,
    target_S0_tau    = target_S0_tau,
    tau              = tau
  )
  
  beta_trt_cal <- calibrate_beta_treatment(
    lp0              = lp0,
    lambda0          = control_cal$lambda0,
    shape            = control_cal$shape,
    target_median_z1 = target_median_z1
  )
  
  lp1    <- lp0 + beta_trt_cal
  S0_tau <- weibull_survival_marginal(tau, lp0, control_cal$lambda0, control_cal$shape)
  S1_tau <- weibull_survival_marginal(tau, lp1, control_cal$lambda0, control_cal$shape)
  
  list(
    lambda0           = control_cal$lambda0,
    weibull_shape     = control_cal$shape,
    beta_treatment_OS = beta_trt_cal,
    checks = data.frame(
      median_z0_target  = target_median_z0,
      median_z0_check   = weibull_median_marginal(lp0, control_cal$lambda0, control_cal$shape),
      S0_tau_target     = target_S0_tau,
      S0_tau_check      = S0_tau,
      median_z1_target  = target_median_z1,
      median_z1_check   = weibull_median_marginal(lp1, control_cal$lambda0, control_cal$shape),
      S1_tau_check      = S1_tau,
      deltaS_tau_check  = S1_tau - S0_tau,
      HR_conditional    = exp(beta_trt_cal),
      logHR_conditional = beta_trt_cal
    )
  )
}


# -----------------------------------------------------------------------------
# Superpopulation
# -----------------------------------------------------------------------------

# Generate a superpopulation with potential outcomes under both arms
#
# Each individual receives latent event times under control (T_true_z0) and
# under treatment (T_true_z1) from independent uniforms, and a common
# administrative censoring time C ~ U(C_low, C_upp) applied to both.
#
# @param N Number of individuals.
# @param C_low,C_upp Bounds of the uniform censoring distribution (months).
# @param lambda0,weibull_shape,beta_treatment_OS Calibrated values, from
#   calibrate_weibull_DGP().
# @return data.frame with covariates, linear predictors, latent and observed
#   times for both arms; `time` and `event` default to the control arm.
generate_superpopulation <- function(
    N,
    C_low,
    C_upp,
    lambda0,
    weibull_shape,
    beta_treatment_OS,
    p_M,
    p_B,
    p_sexF,
    p_resection_0_1_2_3,
    p_trt_baseline_0_1_2,
    p_WHOm2,
    p_loc_superficial_0,
    p_loc_superficial_1,
    p_loc_superficial_2,
    p_loc_superficial_3,
    p_loc_superficial_4,
    p_loc_superficial_5,
    p_multifocal,
    mean_tumor_size,
    sd_tumor_size,
    p_pTERT,
    p_EGFR,
    p_g7l10,
    mean_age,
    sd_age,
    Unmeas_mean,
    Unmeas_sd,
    p_unnec,
    p_unmeas,
    rho_B_TS,
    rho_U_TS,
    int_OS,
    beta_M_OS,
    beta_B_OS,
    beta_age_OS,
    beta_WHOm2_OS,
    beta_sex_OS,
    beta_resection_OS,
    beta_trt_baseline_OS,
    beta_localization_OS,
    beta_multifocal_OS,
    beta_pTERT_OS,
    beta_tumor_size_OS,
    beta_EGFR_OS,
    beta_g7l10_OS,
    beta_Unmeas_OS,
    beta_unmeas_binary
) {
  pop <- generate_covariates(
    N                    = N,
    p_M                  = p_M,
    p_B                  = p_B,
    p_sexF               = p_sexF,
    p_resection_0_1_2_3  = p_resection_0_1_2_3,
    p_trt_baseline_0_1_2 = p_trt_baseline_0_1_2,
    p_WHOm2              = p_WHOm2,
    p_loc_superficial_0  = p_loc_superficial_0,
    p_loc_superficial_1  = p_loc_superficial_1,
    p_loc_superficial_2  = p_loc_superficial_2,
    p_loc_superficial_3  = p_loc_superficial_3,
    p_loc_superficial_4  = p_loc_superficial_4,
    p_loc_superficial_5  = p_loc_superficial_5,
    p_multifocal         = p_multifocal,
    mean_tumor_size      = mean_tumor_size,
    sd_tumor_size        = sd_tumor_size,
    p_pTERT              = p_pTERT,
    p_EGFR               = p_EGFR,
    p_g7l10              = p_g7l10,
    mean_age             = mean_age,
    sd_age               = sd_age,
    Unmeas_mean          = Unmeas_mean,
    Unmeas_sd            = Unmeas_sd,
    p_unnec              = p_unnec,
    p_unmeas             = p_unmeas,
    rho_B_TS             = rho_B_TS,
    rho_U_TS             = rho_U_TS
  )
  
  pop$lp_OS_z0 <- compute_lp_OS(
    dat                  = pop,
    int_OS               = int_OS,
    beta_M_OS            = beta_M_OS,
    beta_B_OS            = beta_B_OS,
    beta_age_OS          = beta_age_OS,
    beta_WHOm2_OS        = beta_WHOm2_OS,
    beta_sex_OS          = beta_sex_OS,
    beta_resection_OS    = beta_resection_OS,
    beta_trt_baseline_OS = beta_trt_baseline_OS,
    beta_localization_OS = beta_localization_OS,
    beta_multifocal_OS   = beta_multifocal_OS,
    beta_pTERT_OS        = beta_pTERT_OS,
    beta_tumor_size_OS   = beta_tumor_size_OS,
    beta_EGFR_OS         = beta_EGFR_OS,
    beta_g7l10_OS        = beta_g7l10_OS,
    beta_Unmeas_OS       = beta_Unmeas_OS,
    beta_unmeas_binary   = beta_unmeas_binary
  )
  pop$lp_OS_z1 <- pop$lp_OS_z0 + beta_treatment_OS
  
  U0 <- runif(N)
  U1 <- runif(N)
  
  pop$T_true_z0 <- (-log(U0) / (lambda0 * exp(pop$lp_OS_z0)))^(1 / weibull_shape)
  pop$T_true_z1 <- (-log(U1) / (lambda0 * exp(pop$lp_OS_z1)))^(1 / weibull_shape)
  
  C <- runif(N, C_low, C_upp)
  
  pop$time_z0  <- pmin(pop$T_true_z0, C)
  pop$event_z0 <- as.numeric(pop$T_true_z0 <= C)
  pop$time_z1  <- pmin(pop$T_true_z1, C)
  pop$event_z1 <- as.numeric(pop$T_true_z1 <= C)
  
  pop$time          <- pop$time_z0
  pop$event         <- pop$event_z0
  pop$resection     <- as.factor(pop$resection)
  pop$trt_type      <- as.factor(pop$trt_type)
  pop$unmeas_binary <- NULL
  
  pop
}


#-----------------------------------------------------------------------------
  # Argument list
  # -----------------------------------------------------------------------------

# Build the list of DGP arguments from one row of workbook parameters
#
# Single place where workbook columns are mapped to function arguments.
# Use with do.call(), e.g.
#   do.call(generate_superpopulation, c(list(N = ..., ...), make_dgp_args(pars)))
#
# @param pars One-row data frame, already checked with check_dgp_params().
# @return Named list of arguments shared by calibrate_weibull_DGP() and
#   generate_superpopulation().
make_dgp_args <- function(pars) {
  list(
    p_M                  = pars$p_M,
    p_B                  = pars$p_B,
    p_sexF               = pars$p_sexF,
    p_resection_0_1_2_3  = c(pars$p_resection_0, pars$p_resection_1,
                             pars$p_resection_2, pars$p_resection_3),
    p_trt_baseline_0_1_2 = c(pars$p_trt_baseline_0, pars$p_trt_baseline_1,
                             pars$p_trt_baseline_2),
    p_WHOm2              = pars$p_WHOm2,
    p_loc_superficial_0  = pars$p_loc_superficial_0,
    p_loc_superficial_1  = pars$p_loc_superficial_1,
    p_loc_superficial_2  = pars$p_loc_superficial_2,
    p_loc_superficial_3  = pars$p_loc_superficial_3,
    p_loc_superficial_4  = pars$p_loc_superficial_4,
    p_loc_superficial_5  = pars$p_loc_superficial_5,
    p_multifocal         = pars$p_multifocal,
    mean_tumor_size      = pars$mean_tumor_size,
    sd_tumor_size        = pars$sd_tumor_size,
    p_pTERT              = pars$p_pTERT,
    p_EGFR               = pars$p_EGFR,
    p_g7l10              = pars$p_g7l10,
    mean_age             = pars$mean_age,
    sd_age               = pars$sd_age,
    Unmeas_mean          = pars$Unmeas_mean,
    Unmeas_sd            = pars$Unmeas_sd,
    p_unnec              = pars$p_unnec,
    p_unmeas             = pars$p_unmeas,
    rho_B_TS             = pars$rho_B_TS,
    rho_U_TS             = pars$rho_U_TS,
    int_OS               = pars$int_OS,
    beta_M_OS            = pars$beta_M_OS,
    beta_B_OS            = pars$beta_B_OS,
    beta_age_OS          = pars$beta_age_OS,
    beta_WHOm2_OS        = pars$beta_WHOm2_OS,
    beta_sex_OS          = pars$beta_sex_OS,
    beta_resection_OS    = pars$beta_resection_OS,
    beta_trt_baseline_OS = pars$beta_trt_baseline_OS,
    beta_localization_OS = pars$beta_localization_OS,
    beta_multifocal_OS   = pars$beta_multifocal_OS,
    beta_pTERT_OS        = pars$beta_pTERT_OS,
    beta_tumor_size_OS   = pars$beta_tumor_size_OS,
    beta_EGFR_OS         = pars$beta_EGFR_OS,
    beta_g7l10_OS        = pars$beta_g7l10_OS,
    beta_Unmeas_OS       = pars$beta_Unmeas_OS,
    beta_unmeas_binary   = pars$beta_unmeas_binary
  )
}


compute_truth_superpop <- function(registry, trial_target, covs, tau) {
  
  reg0   <- registry     |> mutate(Z = 0, S = 0, time = time_z0, event = event_z0)
  trial1 <- trial_target |> mutate(Z = 1, S = 1, time = time_z1, event = event_z1)
  data_full <- bind_rows(reg0, trial1)
  
  S0_true_tau <- mean(trial_target$T_true_z0 > tau, na.rm = TRUE)
  S1_true_tau <- mean(trial_target$T_true_z1 > tau, na.rm = TRUE)
  RMST0_true  <- mean(pmin(trial_target$T_true_z0, tau), na.rm = TRUE)
  RMST1_true  <- mean(pmin(trial_target$T_true_z1, tau), na.rm = TRUE)
  
  # Unadjusted
  fit_unadj <- coxph(Surv(time, event) ~ Z, data = data_full)
  
  # Prognostic score
  fit_prog <- coxph(
    as.formula(paste("Surv(time_z0, event_z0) ~", paste(covs, collapse = " + "))),
    data = registry
  )
  data_full$prog_score <- predict(fit_prog, newdata = data_full, type = "lp")
  fit_adj_prog <- coxph(Surv(time, event) ~ Z + prog_score, data = data_full)
  
  # Propensity score
  fit_ps <- glm(as.formula(paste("S ~", paste(covs, collapse = " + "))),
                data = data_full, family = binomial())
  data_full$prop_score <- predict(fit_ps, newdata = data_full, type = "response")
  fit_adj_prop <- coxph(Surv(time, event) ~ Z + prop_score, data = data_full)
  
  tibble(
    beta_unadj_star    = unname(coef(fit_unadj)["Z"]),
    beta_adj_prog_star = unname(coef(fit_adj_prog)["Z"]),
    beta_adj_prop_star = unname(coef(fit_adj_prop)["Z"]),
    S0_true_tau        = S0_true_tau,
    S1_true_tau        = S1_true_tau,
    deltaS_true        = S1_true_tau - S0_true_tau,
    RMST0_true         = RMST0_true,
    RMST1_true         = RMST1_true,
    deltaRMST_true     = RMST1_true - RMST0_true,
    median_surv0_tau   = median(trial_target$T_true_z0, na.rm = TRUE),
    median_surv1_tau   = median(trial_target$T_true_z1, na.rm = TRUE)
  )
}


