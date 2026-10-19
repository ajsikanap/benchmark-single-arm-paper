# =============================================================================
#
# Selection model for scenario S4 (limited overlap): the trial pool is drawn
# from the superpopulation with probability proportional to exp(lp_sel), where
# lp_sel is a linear selection score whose coefficients (log-odds of trial vs
# registry membership) are defined in config/parameters.R
# (S4_SELECTION_COEFS). Positive coefficients increase the probability of
# entering the trial.
#
# The coefficients were chosen from the target prevalences in the S4 sheet of
# config/simulation_structure.xlsx, attenuated to keep a clinically plausible
# overlap between registry and trial (see Supplementary Table A2).
# =============================================================================


# Design columns of the S4 selection score
#
# Columns are returned in the order of S4_SELECTION_COEFS. Age enters
# standardised on the population passed in.
#
#  dat Superpopulation (resection and trt_type may be factors).
# # Numeric matrix, one column per selection term.
s4_selection_design <- function(dat) {
  resection <- as.integer(as.character(dat$resection))
  trt_type  <- as.integer(as.character(dat$trt_type))
  
  cbind(
    age_z          = scale(dat$age)[, 1],
    WHOm2          = dat$WHOm2,
    resection_3    = as.numeric(resection == 3),
    resection_1    = as.numeric(resection == 1),
    resection_0    = as.numeric(resection == 0),
    trt_type_2     = as.numeric(trt_type == 2),
    trt_type_0     = as.numeric(trt_type == 0),
    trt_type_1     = as.numeric(trt_type == 1),
    multifocal     = dat$multifocal,
    localization_0 = dat$localization_0,
    localization_1 = dat$localization_1,
    localization_2 = dat$localization_2,
    localization_3 = dat$localization_3,
    localization_4 = dat$localization_4,
    pTERT          = dat$pTERT,
    EGFR           = dat$EGFR
  )
}


# Generate the S4 trial pool
#
# Draws N_trial_pool individuals without replacement from the
# superpopulation, with weights exp(lp_sel). Patients with cerebellar
# localization (localization_4 = 1) are excluded, reflecting the trial
# eligibility criteria.
#
# NOTE: the pool is drawn from the whole superpopulation, including the
# individuals later assigned to the registry (02_generate_populations.R), so
# the S4 trial pool and the registry can share individuals.
#
#  superpop     Superpopulation from generate_superpopulation().
#  N_trial_pool Size of the trial pool.
#  coefs        Named vector of selection coefficients
#                     (S4_SELECTION_COEFS), in the order of s4_selection_design().
#  seed         Seed set before sampling.
# # Trial pool, with the selection score in column s4_selection_lp and
#   resection / trt_type as integers.
make_s4_trial_pool <- function(superpop, N_trial_pool, coefs, seed) {
  X <- s4_selection_design(superpop)
  stopifnot(identical(colnames(X), names(coefs)))
  
  set.seed(seed)
  
  # Terms are accumulated one at a time, in a fixed order, so that the score
  # is reproducible to the last digit.
  lp <- coefs[[1]] * X[, 1]
  for (j in seq_along(coefs)[-1]) lp <- lp + coefs[[j]] * X[, j]
  
  # resection and trt_type are returned as integers, as in the original
  # pipeline; they are converted back to factors when the scenario pools are
  # built.
  superpop$resection       <- as.integer(as.character(superpop$resection))
  superpop$trt_type        <- as.integer(as.character(superpop$trt_type))
  superpop$s4_selection_lp <- lp
  
  superpop |>
    dplyr::filter(localization_4 == 0) |>
    dplyr::slice_sample(
      n         = N_trial_pool,
      weight_by = exp(s4_selection_lp),
      replace   = FALSE
    )
}