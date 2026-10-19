# =============================================================================
# R/tables_dgp.R
#
# Tables describing the data-generating process:
#   - table_dgp_parameters():     covariate distributions and prognostic
#                                 coefficients (Table 1)
#   - table_weibull_calibration(): calibrated Weibull parameters and checks
#                                 (Supplementary Table A1)
#   - table_s4_shift():           S4 covariate shift, target vs observed
#                                 (Supplementary Table A2)
#
# All functions return numeric tibbles (NA where a cell does not apply);
# formatting is left to the calling document.
# =============================================================================


# Labels and grouping shared by the tables
dgp_covariates <- tibble::tribble(
  ~key,             ~section,                 ~label,
  "B",              "Molecular markers",      "Brain penetrance (B)",
  "M",              "Molecular markers",      "MGMT methylation (M)",
  "pTERT",          "Molecular markers",      "pTERT mutation",
  "EGFR",           "Molecular markers",      "EGFR alteration",
  "g7l10",          "Molecular markers",      "Chr7+/Chr10\u2212",
  "age",            "Demographics",           "Age (years)",
  "sex",            "Demographics",           "Sex (female)",
  "WHOm2",          "Demographics",           "WHO PS \u2265 2 at diagnosis",
  "resection_0",    "Extent of resection",    "0: Biopsy",
  "resection_1",    "Extent of resection",    "1: Partial resection",
  "resection_2",    "Extent of resection",    "2: Subtotal resection",
  "resection_3",    "Extent of resection",    "3: Total (>95%)",
  "trt_type_0",     "Treatment intensity",    "0: Surgery only",
  "trt_type_1",     "Treatment intensity",    "1: Intermediate",
  "trt_type_2",     "Treatment intensity",    "2: Aggressive",
  "tumor_size",     "Tumour characteristics", "Tumour size (mL)",
  "multifocal",     "Tumour characteristics", "Multifocal disease",
  "localization_0", "Tumour localization",    "Frontal",
  "localization_1", "Tumour localization",    "Parietal",
  "localization_2", "Tumour localization",    "Temporal",
  "localization_3", "Tumour localization",    "Occipital",
  "localization_4", "Tumour localization",    "Cerebellum",
  "localization_5", "Tumour localization",    "Brainstem"
)


# Target prevalence (or mean / SD) of each covariate in one row of parameters
target_summary <- function(pars) {
  tibble::tibble(
    key   = dgp_covariates$key,
    value = c(
      pars$p_B, pars$p_M, pars$p_pTERT, pars$p_EGFR, pars$p_g7l10,
      pars$mean_age, pars$p_sexF, pars$p_WHOm2,
      pars$p_resection_0, pars$p_resection_1, pars$p_resection_2, pars$p_resection_3,
      pars$p_trt_baseline_0, pars$p_trt_baseline_1, pars$p_trt_baseline_2,
      pars$mean_tumor_size, pars$p_multifocal,
      pars$p_loc_superficial_0, pars$p_loc_superficial_1, pars$p_loc_superficial_2,
      pars$p_loc_superficial_3, pars$p_loc_superficial_4, pars$p_loc_superficial_5
    ),
    sd = ifelse(key == "age", pars$sd_age,
                ifelse(key == "tumor_size", pars$sd_tumor_size, NA_real_))
  )
}


# Observed prevalence (or mean / SD) of each covariate in a population
observed_summary <- function(dat) {
  res <- as.integer(as.character(dat$resection))
  trt <- as.integer(as.character(dat$trt_type))
  tibble::tibble(
    key   = dgp_covariates$key,
    value = c(
      mean(dat$B), mean(dat$M), mean(dat$pTERT), mean(dat$EGFR), mean(dat$g7l10),
      mean(dat$age), mean(dat$sex), mean(dat$WHOm2),
      mean(res == 0), mean(res == 1), mean(res == 2), mean(res == 3),
      mean(trt == 0), mean(trt == 1), mean(trt == 2),
      mean(dat$tumor_size), mean(dat$multifocal),
      mean(dat$localization_0), mean(dat$localization_1), mean(dat$localization_2),
      mean(dat$localization_3), mean(dat$localization_4), mean(dat$localization_5)
    ),
    sd = ifelse(key == "age", sd(dat$age),
                ifelse(key == "tumor_size", sd(dat$tumor_size), NA_real_))
  )
}


# Table 1: covariate distributions and prognostic coefficients
#
# Extent of resection and treatment intensity are ordinal with a linear
# effect per level: the log-HR is reported once per variable (per level),
# on the row of the highest level the contrast with the reference
# (3 x per-level log-HR for resection, 2 x for treatment).
#
# pars       One row of registry parameters (check_dgp_params()).
# references Named character vector of references, by covariate key.
#  tibble: section, label, value, sd, log_HR, HR, note, reference.
table_dgp_parameters <- function(pars, references) {
  beta <- c(
    B = pars$beta_B_OS, M = pars$beta_M_OS, pTERT = pars$beta_pTERT_OS,
    EGFR = pars$beta_EGFR_OS, g7l10 = pars$beta_g7l10_OS,
    age = pars$beta_age_OS, sex = pars$beta_sex_OS, WHOm2 = pars$beta_WHOm2_OS,
    resection_0 = 0,
    resection_1 = 1 * pars$beta_resection_OS,
    resection_2 = 2 * pars$beta_resection_OS,
    resection_3 = 3 * pars$beta_resection_OS,
    trt_type_0 = 0,
    trt_type_1 = 1 * pars$beta_trt_baseline_OS,
    trt_type_2 = 2 * pars$beta_trt_baseline_OS,
    tumor_size = pars$beta_tumor_size_OS, multifocal = pars$beta_multifocal_OS,
    localization_0 = pars$beta_localization_OS, localization_1 = pars$beta_localization_OS,
    localization_2 = pars$beta_localization_OS, localization_3 = pars$beta_localization_OS,
    localization_4 = pars$beta_localization_OS, localization_5 = pars$beta_localization_OS
  )
  
  note <- c(
    age        = "per year, centred at 60",
    tumor_size = "per mL",
    resection_0 = "reference",
    resection_1 = "linear per level",
    resection_2 = "linear per level",
    resection_3 = "linear per level",
    trt_type_0 = "reference",
    trt_type_1 = "linear per level",
    trt_type_2 = "linear per level"
  )
  
  dgp_covariates |>
    dplyr::left_join(target_summary(pars), by = "key") |>
    dplyr::mutate(
      log_HR    = unname(beta[key]),
      HR        = exp(log_HR),
      note      = unname(note[key]),
      reference = unname(references[key])
    ) |>
    dplyr::select(section, label, value, sd, log_HR, HR, note, reference)
}


# Supplementary Table A1: calibrated Weibull parameters and checks
#
# Calibrated values are the marginal quantities on the calibration
# population; verified values are empirical quantities of the latent event
# times in the generated superpopulation.
#
# cal      Output of calibrate_weibull_DGP().
# superpop Output of generate_superpopulation().
# tau      Landmark time (months).
#  tibble: parameter, target, calibrated, verified.
table_weibull_calibration <- function(cal, superpop, tau) {
  ch <- cal$checks
  S0_emp <- mean(superpop$T_true_z0 > tau)
  S1_emp <- mean(superpop$T_true_z1 > tau)
  
  tibble::tibble(
    parameter = c(
      "Weibull scale (\u03bb\u2080)",
      "Weibull shape (k)",
      "Control median (months)",
      sprintf("Control S(%g)", tau),
      "Treatment median (months)",
      sprintf("Treatment S(%g)", tau),
      sprintf("\u0394S(%g)", tau),
      "Treatment log-HR (\u03b2_Z, conditional)",
      "Treatment HR (conditional)"
    ),
    target = c(
      NA, NA, ch$median_z0_target, ch$S0_tau_target, ch$median_z1_target,
      NA, NA, NA, NA
    ),
    calibrated = c(
      cal$lambda0, cal$weibull_shape, ch$median_z0_check, ch$S0_tau_check,
      ch$median_z1_check, ch$S1_tau_check, ch$deltaS_tau_check,
      cal$beta_treatment_OS, exp(cal$beta_treatment_OS)
    ),
    verified = c(
      NA, NA,
      stats::median(superpop$T_true_z0), S0_emp,
      stats::median(superpop$T_true_z1), S1_emp,
      S1_emp - S0_emp, NA, NA
    )
  )
}


# Supplementary Table A2: S4 covariate shift, target vs observed
#
# reg_target   Registry target parameters (S4 sheet, registry row).
# trial_target Trial target parameters (S4 sheet, trial row).
# superpop     Superpopulation (observed registry distribution).
# trial_pool   S4 trial pool, before restriction to B = 1.
# coefs        S4 selection coefficients (S4_SELECTION_COEFS).
# references   Named character vector of references, by covariate key.
#  tibble: section, label, registry/trial target and observed values
#   (mean for continuous covariates), selection log-OR, note, reference.
table_s4_shift <- function(reg_target, trial_target, superpop, trial_pool,
                           coefs, references) {
  sel <- c(coefs, resection_2 = 0)
  sel_key <- c(age = "age_z")
  
  note <- c(
    B = "not in the selection model", M = "not in the selection model",
    g7l10 = "not in the selection model", sex = "not in the selection model",
    tumor_size = "not in the selection model",
    localization_5 = "not in the selection model",
    age = "per SD of age",
    resection_2 = "reference",
    localization_4 = "excluded from the trial pool"
  )
  
  dgp_covariates |>
    dplyr::mutate(
      registry_target = target_summary(reg_target)$value,
      trial_target    = target_summary(trial_target)$value,
      registry_obs    = observed_summary(superpop)$value,
      trial_obs       = observed_summary(trial_pool)$value,
      sel_name        = ifelse(key %in% names(sel_key), sel_key[key], key),
      selection_logOR = ifelse(sel_name %in% names(sel), unname(sel[sel_name]), 0),
      note            = unname(note[key]),
      reference       = unname(references[key])
    ) |>
    dplyr::select(section, label, registry_target, trial_target,
                  registry_obs, trial_obs, selection_logOR, note, reference)
}