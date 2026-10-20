# =============================================================================
# 04_run_simulation.R
#
# Main simulation study, scenarios S0-S5.
# Produces ONLY the raw per-iteration results; every performance measure
# (bias, RMSE, type I error, power, MCSE) is computed downstream in
# 06_paper_outputs.Rmd.
#
# Design notes
# - Parallelised over scenarios (one scenario per worker) with future/furrr.
# - Random numbers: iteration i uses seeds[i] in every scenario and for every
#   registry size (registries are nested within an iteration), so comparisons
#   across scenarios and across N are paired (common random numbers).
#   Seeding is handled inside run_scenario_analysis()/one_resample2();
#   furrr's own seeding is therefore disabled (seed = NULL).
# - Each run returns both "null" and "alternative" rows (column
#   effect_setting): one_resample2v2 always evaluates both, so the
#   effect_setting argument of run_scenario_analysis() is not used here.
# - run_scenario_analysis() reads `seeds` and `trial_pool_scenarios` from the
#   global environment; furrr exports them to the workers automatically.
# - Each scenario is checkpointed to results/checkpoints_<N_BOOT>/ together
#   with the parameters that generated it; a checkpoint is reused only if
#   those parameters are identical to the current ones.
# =============================================================================

library(survival)
library(MatchIt)
library(dplyr)
library(purrr)
library(future)
library(furrr)


# ── Functions and inputs ──────────────────────────────────────────────────────
source(here::here("functions/load_functions.R"))
source(here::here("config/parameters.R"))


# ── Global parameters ─────────────────────────────────────────────────────────
N_BOOT         <- 10
N_REGISTRY     <- N_REG_VALUES

scenario1 <- readxl::read_excel(here::here("config/simulation_structure.xlsx"), sheet = "S1")
scenario4 <- readxl::read_excel(here::here("config/simulation_structure.xlsx"), sheet = "S4")

seeds <- read.table(here::here("config/seeds.txt"), header = FALSE)$V1

load(here::here("data/registry_generated.rda"))                  # registry_generated
load(here::here("data/trial_pool_scenarios.rda"))                # trial_pool_scenarios
load(here::here("data/truth_scenarios_superpopulation.rda"))     # truth_scenarios

# ── Scenario configurations ───────────────────────────────────────────────────
# Only the arguments that differ across scenarios are spelled out.
make_cfg <- function(scenario, label, sc,
                     B_setting = "all", include_B_in_score = FALSE) {
  list(
    label                         = label,
    scenario_name                 = scenario,
    trial_pool_name               = scenario,
    sc                            = sc,
    B_setting                     = B_setting,
    include_B_in_score            = include_B_in_score,
    include_unneccessary_in_score = FALSE
  )
}

scenario_configs <- list(
  S0 = make_cfg("S0", "S0_unfiltered_oracle",     scenario1, B_setting = "unfiltered", include_B_in_score = TRUE),
  S1 = make_cfg("S1", "S1_B_known_registry_all",  scenario1, include_B_in_score = TRUE),
  S2 = make_cfg("S2", "S2_B_unknown",             scenario1),
  S3A = make_cfg("S3A", "S3A_B_unknown_M_shift",     scenario1),
  S3B = make_cfg("S3B", "S3B_unfavorable_cov_shift", scenario1),
  S4 = make_cfg("S4", "S4_limited_overlap",       scenario4)
)

# ── Pre-run checks ────────────────────────────────────────────────────────────
stopifnot(
  "seeds.txt has fewer than N_BOOT seeds" = length(seeds) >= N_BOOT,
  "seeds.txt contains NA"                 = !anyNA(seeds[seq_len(N_BOOT)])
)

for (sc_name in names(scenario_configs)) {
  cfg <- scenario_configs[[sc_name]]
  if (!cfg$trial_pool_name %in% names(trial_pool_scenarios)) {
    stop(sprintf("Trial pool '%s' not found in trial_pool_scenarios (scenario %s).",
                 cfg$trial_pool_name, sc_name))
  }
  if (sum(truth_scenarios$scenario == cfg$scenario_name) != 1) {
    stop(sprintf("Scenario '%s' must have exactly one row in truth_scenarios.",
                 cfg$scenario_name))
  }
}

# ── Checkpoints ───────────────────────────────────────────────────────────────
run_params <- list(N_BOOT = N_BOOT, N_TRIAL = N_TRIAL, TAU = TAU,
                   RATIO_CONTROLS = RATIO_CONTROLS, N_REGISTRY = N_REGISTRY)

ckpt_dir <- here::here("results", sprintf("checkpoints_%d", N_BOOT))
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

ckpt_file <- function(sc_name) file.path(ckpt_dir, paste0(sc_name, ".rds"))

run_one_scenario <- function(sc_name, cfg) {
  file        <- ckpt_file(sc_name)
  fingerprint <- list(params = run_params, cfg = cfg)
  
  if (file.exists(file)) {
    ck <- readRDS(file)
    if (is.list(ck) && identical(ck$fingerprint, fingerprint)) {
      message(sprintf("[cached] %s", sc_name))
      return(data.frame(scenario = sc_name, status = "cached", error_msg = NA_character_))
    }
    message(sprintf("[rerun] %s: checkpoint built with different parameters", sc_name))
  }
  
  tryCatch({
    res <- run_scenario_analysis(
      label                         = cfg$label,
      scenario_name                 = cfg$scenario_name,
      trial_pool_name               = cfg$trial_pool_name,
      sc                            = cfg$sc,
      B_setting                     = cfg$B_setting,
      include_B_in_score            = cfg$include_B_in_score,
      include_unneccessary_in_score = cfg$include_unneccessary_in_score,
      registry_pool                 = registry_generated,
      truth_scenarios               = truth_scenarios,
      n_reg_values                  = N_REGISTRY,
      n_boot                        = N_BOOT,
      n_trial                       = N_TRIAL,
      ratio_controls                = RATIO_CONTROLS,
      tau                           = TAU
    )
    # Atomic write: a crash during saveRDS cannot leave a truncated checkpoint.
    tmp <- paste0(file, ".tmp")
    saveRDS(list(fingerprint = fingerprint, res = res), tmp)
    file.rename(tmp, file)
    message(sprintf("[ok] %s", sc_name))
    data.frame(scenario = sc_name, status = "ok", error_msg = NA_character_)
  },
  error = function(e) {
    message(sprintf("[ERROR] %s: %s", sc_name, conditionMessage(e)))
    data.frame(scenario = sc_name, status = "error", error_msg = conditionMessage(e))
  })
}

# ── Run ───────────────────────────────────────────────────────────────────────
n_workers <- min(length(scenario_configs), max(1L, parallelly::availableCores() - 1L))
plan(multisession, workers = n_workers)
message(sprintf("Running %d scenarios on %d workers.", length(scenario_configs), n_workers))

run_log <- furrr::future_imap(
  scenario_configs,
  \(cfg, sc_name) run_one_scenario(sc_name, cfg),
  .options = furrr::furrr_options(seed = NULL)
) |> dplyr::bind_rows()

plan(sequential)
print(run_log)

if (any(run_log$status == "error")) {
  stop("Scenario(s) failed: ",
       paste(run_log$scenario[run_log$status == "error"], collapse = ", "),
       ". Nothing saved.")
}

# ── Collect results from the checkpoints of the current scenarios ────────────
res_boot <- purrr::map(names(scenario_configs), \(s) readRDS(ckpt_file(s))$res) |>
  dplyr::bind_rows()

# Iterations that crashed (status == "error") carry no effect_setting and would
# be dropped silently by the null/alternative split: keep them apart and report.
# These are code failures, distinct from match_failed (degenerate matching),
# which stays in res_boot as a reportable outcome.
res_boot_errors <- res_boot |> dplyr::filter(status == "error")
res_boot        <- res_boot |> dplyr::filter(status == "ok")

if (nrow(res_boot_errors) > 0) {
  warning(sprintf("%d iteration(s) crashed; see res_boot_errors.", nrow(res_boot_errors)))
  print(dplyr::count(res_boot_errors, scenario, registry_n, error_msg))
}

stopifnot(
  "effect_setting missing in some rows" =
    all(res_boot$effect_setting %in% c("null", "alternative"))
)

# Every scenario x registry size x effect setting should contain N_BOOT iterations.
iter_check <- res_boot |>
  dplyr::distinct(scenario, registry_n, effect_setting, iter) |>
  dplyr::count(scenario, registry_n, effect_setting, name = "n_iter") |>
  dplyr::filter(n_iter != N_BOOT)
if (nrow(iter_check) > 0) {
  warning("Some cells have fewer than N_BOOT iterations:")
  print(iter_check)
}

res_boot_alt  <- res_boot |> dplyr::filter(effect_setting == "alternative")
res_boot_null <- res_boot |> dplyr::filter(effect_setting == "null")

# ── Save ──────────────────────────────────────────────────────────────────────
output_file <- here::here(sprintf("results/results_raw_%d.RData", N_BOOT))

save(
  res_boot_alt, res_boot_null, res_boot_errors,
  truth_scenarios, scenario_configs,
  N_BOOT, N_TRIAL, TAU, RATIO_CONTROLS, N_REGISTRY, ALPHA,
  file = output_file
)

message(sprintf("Done. Results saved to %s", output_file))
