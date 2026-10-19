# =============================================================================
# analysis/main/02_generate_populations.R
#
# Builds the registry and the scenario-specific trial pools (S0-S4) from the
# superpopulation, and computes the true estimands in each trial pool.
#
# Input : data/superpopulation_generated_weibull.rda   (01_dgp_weibull.Rmd)
#         data/trial_pool_s4_generated_weibull_s4.rda  (01_dgp_weibull.Rmd)
# Output: data/registry_generated.rda              registry (S = 0)
#         data/trial_pool_scenarios.rda            list of trial pools, by scenario
#         data/truth_scenarios_superpopulation.rda true estimands, by scenario
# Time  : about 1 minute
#
# As in Table 3: Scenarios (trial pools of N_TRIAL_POOL individuals). Internal codes are
# used in data and results; paper labels in brackets:
#   S0        no restriction on B
#   S1        B = 1 (B observed in the registry)
#   S2        B = 1 (B unobserved: same pool as S1, differs only in the analysis)
#   S3 (S3A)  B = 1 and favourable shift in M
#   S4        B = 1 and selection model on all covariates (limited overlap)
#   S5 (S3B)  B = 1 and unfavourable shift in M, WHO PS and multifocality
#
# Conventions: registry has S = 0 and time / event = control outcome (z0);
# trial pools have S = 1 and time / event = treated outcome (z1), except S4
# (see note below).
# =============================================================================

library(here)
library(dplyr)
library(survival)

source(here("config/parameters.R"))

load(here("data/superpopulation_generated_weibull.rda"))
load(here("data/trial_pool_s4_generated_weibull_s4.rda"))

set.seed(SEED_MASTER)

# The order of the random operations below must not be changed: it
# determines the reproducibility of the pools used in the paper.


# -----------------------------------------------------------------------------
# Registry and trial-eligible population (disjoint)
# -----------------------------------------------------------------------------
ind <- sample(c(TRUE, FALSE), nrow(superpopulation_generated), replace = TRUE,
              prob = c(P_REGISTRY_SPLIT, 1 - P_REGISTRY_SPLIT))

registry_generated <- superpopulation_generated[ind, ] |>
  mutate(time = time_z0, event = event_z0, S = 0)

trial_pool_generated <- superpopulation_generated[!ind, ]


# -----------------------------------------------------------------------------
# S0: no restriction on B
# -----------------------------------------------------------------------------
S0 <- trial_pool_generated |>
  mutate(time = time_z1, event = event_z1, S = 1) |>
  slice_sample(n = N_TRIAL_POOL, replace = FALSE)


# -----------------------------------------------------------------------------
# S1 / S2: restriction to B = 1
# -----------------------------------------------------------------------------
S1 <- trial_pool_generated |>
  filter(B == 1) |>
  mutate(time = time_z1, event = event_z1, S = 1) |>
  slice_sample(n = N_TRIAL_POOL, replace = FALSE)

S2 <- S1


# -----------------------------------------------------------------------------
# S3 (paper S3A) and S5 (paper S3B): B = 1 and covariate shift
#
# Importance resampling: the trial-eligible population is resampled with
# replacement, with weights p_target / p_ref (x = 1) and
# (1 - p_target) / (1 - p_ref) (x = 0), then restricted to B = 1, and
# N_TRIAL_POOL individuals are drawn without replacement. Because of the
# first step, the pool can contain the same individual more than once.
#
# S5 is generated after S4: the order of the random operations is that of the
# original pipeline.
# -----------------------------------------------------------------------------

# S3 (paper S3A): favourable shift in M (reference prevalence: registry target)
w_S3A <- ifelse(trial_pool_generated$M == 1,
               P_M_TARGET_S3A / P_M_REF_S3A,
               (1 - P_M_TARGET_S3A) / (1 - P_M_REF_S3A))

S3A <- trial_pool_generated |>
  slice_sample(n = nrow(trial_pool_generated), weight_by = w_S3A, replace = TRUE) |>
  filter(B == 1) |>
  mutate(time = time_z1, event = event_z1, S = 1) |>
  slice_sample(n = N_TRIAL_POOL, replace = FALSE)


# -----------------------------------------------------------------------------
# S3B: B = 1 and unfavourable shift in M, WHO PS and multifocality
# Reference prevalences: empirical, in the trial-eligible population.
# -----------------------------------------------------------------------------
p_M_ref     <- mean(trial_pool_generated$M)
p_WHO_ref   <- mean(trial_pool_generated$WHOm2)
p_multi_ref <- mean(trial_pool_generated$multifocal)

w_M <- ifelse(trial_pool_generated$M == 1,
              P_M_TARGET_S3B / p_M_ref,
              (1 - P_M_TARGET_S3B) / (1 - p_M_ref))

w_WHO <- ifelse(trial_pool_generated$WHOm2 == 1,
                P_WHO_TARGET_S3B / p_WHO_ref,
                (1 - P_WHO_TARGET_S3B) / (1 - p_WHO_ref))

w_multi <- ifelse(trial_pool_generated$multifocal == 1,
                  P_MULTI_TARGET_S3B / p_multi_ref,
                  (1 - P_MULTI_TARGET_S3B) / (1 - p_multi_ref))

w_S3B <- w_M * w_WHO * w_multi

S3B <- trial_pool_generated |>
  slice_sample(n = nrow(trial_pool_generated), weight_by = w_S3B, replace = TRUE) |>
  filter(B == 1) |>
  mutate(time = time_z1, event = event_z1, S = 1) |>
  slice_sample(n = N_TRIAL_POOL, replace = FALSE)


trial_pool_scenarios <- list(S0 = S0, S1 = S1, S2 = S2, S3A = S3A, S3B = S3B, S4 = S4)



# -----------------------------------------------------------------------------
# S4: B = 1 and selection model (pool generated in 01_dgp_weibull.Rmd)
#
# NOTE 1: time / event are NOT set to the treated outcome, unlike the other
# scenarios: they keep the control outcome inherited from the
# superpopulation. This affects the Cox estimates in one_resample2(), which
# use time / event of the trial.
# NOTE 2: the S4 pool was drawn from the whole superpopulation, so it can
# share individuals with the registry.
# Both kept as in the original pipeline.
# -----------------------------------------------------------------------------
S4 <- trial_pool_s4 |>
  mutate(
    resection = as.factor(resection),
    trt_type  = as.factor(trt_type),
    S         = 1
  ) |>
  select(-s4_selection_lp) |>
  filter(B == 1) |>
  slice_sample(n = N_TRIAL_POOL, replace = FALSE)



# -----------------------------------------------------------------------------
# True estimands
#
# Survival, RMST and median: latent event times of the trial pool.
# Cox log-HRs: observed times of the full registry under control and of the
# trial pool under treatment; unadjusted, adjusted for a prognostic score
# (Cox model on the registry) and for a propensity score (logistic model for
# trial membership).
# -----------------------------------------------------------------------------
truth_scenarios <- purrr::imap_dfr(
  trial_pool_scenarios,
  \(pool, name) compute_truth_superpop(
    registry     = registry_generated,
    trial_target = pool,
    covs         = COVS_SCORE_MAIN,
    tau          = TAU
  ) |>
    mutate(scenario = name, .before = 1)
) 

 # Paper version of Table 3: order, labels and column names (display only)
 S0_ref <- truth_scenarios$S0_true_tau[truth_scenarios$scenario == "S0"]
 
truth_scenarios |>
   filter(scenario != "S2") |>                                   # identical to S1
   mutate(
     theta_true = log(S0_true_tau) / log(S0_ref),
     delta_true = log(theta_true),
     scenario   = factor(scenario, levels = c("S0", "S1", "S3A", "S3B", "S4"))
   ) |>
   arrange(scenario) |>
   mutate(scenario = dplyr::recode(as.character(scenario),
                                   S1 = "S1-S2")) |>
   select(
     "Scenario"         = scenario,
     "S\u2080(t)*"      = S0_true_tau,
     "S\u2081(t)*"      = S1_true_tau,
     "\u0394S(t)*"      = deltaS_true,
     "RMST\u2080*"      = RMST0_true,
     "RMST\u2081*"      = RMST1_true,
     "\u0394RMST*"      = deltaRMST_true,
     "Median\u2080"     = median_surv0_tau,
     "Median\u2081"     = median_surv1_tau,
     # "\u03b2 unadj*"     = beta_unadj_star,
     # "\u03b2 adj (prog)*" = beta_adj_prog_star,
     # "\u03b2 adj (prop)*" = beta_adj_prop_star,
     "\u03b8_true"       = theta_true,
     "\u03b4_true"       = delta_true
   ) |>
   mutate(across(where(is.numeric), \(x) round(x, 3))) |> flextable::flextable()

# -----------------------------------------------------------------------------
# Checks: realised prevalences and distinct individuals in each pool
# -----------------------------------------------------------------------------
shift_check <- purrr::imap_dfr(
  c(list(registry = registry_generated), trial_pool_scenarios),
  \(d, name) tibble(
    population   = name,
    n            = nrow(d),
    n_unique_id  = n_distinct(d$id),
    p_B          = mean(d$B),
    p_M          = mean(d$M),
    p_WHOm2      = mean(d$WHOm2),
    p_multifocal = mean(d$multifocal)
  )
)
print(shift_check, width = Inf)
print(truth_scenarios, width = Inf)


# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
save(registry_generated,   file = here("data/registry_generated.rda"))
save(trial_pool_scenarios, file = here("data/trial_pool_scenarios.rda"))
save(truth_scenarios,      file = here("data/truth_scenarios_superpopulation.rda"))