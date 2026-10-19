# =============================================================================
# config/parameters.R
#
# Single source of truth for all design and simulation parameters.
# Source at the top of every script:
#   source(here::here("config/parameters.R"))
#
# DGP covariate parameters (prevalences, log-HRs, correlations, censoring
# bounds) are NOT defined here: they live in config/simulation_structure.xlsx.
# Rule: every number used in the pipeline is defined exactly once, either in
# this file or in the xlsx, never in a script or a function body.
# =============================================================================


# -- Reproducibility ----------------------------------------------------------
SEED_MASTER <- 1111                  # population generation and calibration
SEED_FILE   <- "config/seeds.txt"    # one seed per Monte Carlo replicate,
# used as seeds[iter] in one_resample2()


# -- Clinical targets for the Weibull DGP calibration -------------------------
TAU              <- 12     # landmark time (months) for S(tau) and RMST
TARGET_MEDIAN_Z0 <- 7.5    # median OS under control (months)
TARGET_S0_TAU    <- 0.32   # S_0(tau) under control
TARGET_MEDIAN_Z1 <- 13     # median OS under treatment (months)


# -- Population sizes ---------------------------------------------------------
# Main analysis (S0-S4)
N_CALIB_MAIN    <- 100000  # calibration population
N_SUPERPOP_MAIN <- 100000  # superpopulation, split into registry / trial pool
N_TRIAL_POOL_S4 <- 50000   # dedicated trial pool for limited overlap (S4)

# Sensitivity analysis (beta_B x rho_B_TS grid)
N_CALIB_SENS    <- 50000
N_SUPERPOP_SENS <- 50000


# -- Registry / trial pool construction ---------------------------------------
P_REGISTRY_SPLIT <- 0.4    # share of the superpopulation assigned to the registry
N_TRIAL_POOL     <- 20000  # size of each scenario-specific trial pool


# -- Covariate shift targets --------------------------------------------------
# Shift on M only (S3A)
P_M_TARGET_S3A <- 0.60
P_M_REF_S3A    <- 0.46     # TO CHECK: hardcoded, whereas S3B computes the
# reference prevalence from the trial pool

# Unfavourable shift on M, WHO PS and multifocality (S3B)
P_M_TARGET_S3B     <- 0.30
P_WHO_TARGET_S3B   <- 0.45
P_MULTI_TARGET_S3B <- 0.50



# -- Single-arm trial design --------------------------------------------------
N_TRIAL        <- 24       # trial sample size
RATIO_CONTROLS <- 5        # matched controls per trial patient
ALPHA          <- 0.10     # one-sided significance level
N_REG_VALUES   <- c(75, 100, 125, 175, 250, 500, 1000)  # registry sizes
C_upp          <- 30
C_low          <- 0


# -- Monte Carlo --------------------------------------------------------------
N_BOOT <- 10000            # set to 50 for verification runs


# -- Sensitivity analysis grid ------------------------------------------------
# Calibration grid (all cells calibrated in the sensitivity DGP)
BETA_B_CALIB_GRID   <- c(-1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1.00)
RHO_B_TS_CALIB_GRID <- c(-0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75)

# Subset of cells actually simulated
BETA_B_SENS     <- c(-1, 0, 1)
RHO_B_TS_SENS   <- c(-0.75, -0.5, 0, 0.5, 0.75)
N_REGISTRY_SENS <- 250


# -- Covariates for the true estimands ----------------------------------------
# Main analysis: B included (oracle set used to compute the truth)
COVS_SCORE_MAIN <- c(
  "M", "B", "age", "sex", "WHOm2",
  "factor(resection)", "factor(trt_type)",
  "multifocal", "pTERT", "tumor_size", "EGFR", "g7l10"
)

# Sensitivity analysis: B excluded; must match
# make_covariates(sc, include_B = FALSE, include_unneccessary = FALSE)
COVS_SCORE_SENS <- c(
  "M", "age", "sex", "WHOm2",
  "resection", "trt_type",
  "multifocal", "pTERT", "tumor_size", "EGFR", "g7l10"
)


# -- Paths --------------------------------------------------------------------
PATH_DGP_XLSX <- "config/simulation_structure.xlsx"


# -- S4 selection model (log-odds of trial vs registry membership) ------------
# Order matters: it must match s4_selection_design() in R/dgp_s4.R.
# Age enters standardised; localization_4 (cerebellum) is also excluded from
# the trial pool by design.
S4_SELECTION_COEFS <- c(
  age_z          = -0.25,
  WHOm2          = -1.40,
  resection_3    =  0.20,
  resection_1    = -0.40,
  resection_0    = -0.10,
  trt_type_2     =  0.60,
  trt_type_0     = -0.50,
  trt_type_1     = -0.30,
  multifocal     = -0.40,
  localization_0 = -0.53,
  localization_1 = -0.48,
  localization_2 = -0.20,
  localization_3 = -1.00,
  localization_4 = -0.50,
  pTERT          =  0.63,
  EGFR           = -0.22
)