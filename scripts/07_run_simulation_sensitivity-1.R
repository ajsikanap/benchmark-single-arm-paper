# =============================================================================
# 07_run_simulation_sensitivity.R
#
# Sensitivity analysis S5, grid (beta_B, rho_B_TS). Same structure as
# 04_main_simulations_S0_S4.R: produces ONLY the raw per-iteration results;
# every performance measure is computed downstream in 08_analyse_sensitivity.Rmd.
#
# Input : data/superpop_results_sensitivity.rda  (02_dgp_sensitivity.Rmd)
#         data/truth_sensitivity.rda             (02_dgp_sensitivity.Rmd)
# Output: results/results_sensitivity_raw_<N_BOOT>.RData
#
# Design notes
# - Each cell is analysed exactly as S2 in 04: B_setting = "all", B not in the
#   scores, same functions, same seeds (seeds[i] at iteration i in every cell:
#   common random numbers across cells and with the main analysis).
# - Registry sizes: N_REGISTRY_SENS (reported in the paper) and
#   max(N_REG_VALUES). make_nested_registries() draws max(n) individuals and
#   takes the first n: with the same maximum as in 04, the registries of size
#   N_REGISTRY_SENS coincide with those of the main analysis, so the control
#   cell reproduces S2 iteration by iteration. The largest size also shows
#   whether the plateau of the matching error holds across the grid.
# - Parallelised over cells (one cell per worker) with future/furrr; seeding
#   is handled inside run_scenario_analysis()/one_resample2(), so furrr's own
#   seeding is disabled (seed = NULL).
# - run_scenario_analysis() reads `seeds` and `trial_pool_scenarios` from the
#   global environment: run_one_cell() assigns both explicitly in the global
#   environment of the worker before calling it.
# - Each cell is checkpointed to results/checkpoints_sensitivity_<N_BOOT>/
#   with the parameters that generated it; a checkpoint is reused only if
#   they are identical to the current ones.
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

N_REGISTRY <- sort(unique(c(N_REGISTRY_SENS, max(N_REG_VALUES))))

scenario1 <- readxl::read_excel(here::here(PATH_DGP_XLSX), sheet = "S1")

seeds <- read.table(here::here(SEED_FILE), header = FALSE)$V1

load(here::here("data/superpop_results_sensitivity.rda"))   # superpop_results
load(here::here("data/truth_sensitivity.rda"))              # truth_sensitivity


# ── Pre-run checks ────────────────────────────────────────────────────────────
covs_analysis <- make_covariates(scenario1, include_B = FALSE, include_unneccessary = FALSE)

stopifnot(
  "seeds.txt has fewer than N_BOOT seeds" = length(seeds) >= N_BOOT,
  "seeds.txt contains NA"                 = !anyNA(seeds[seq_len(N_BOOT)]),
  "analysis covariates differ from COVS_SCORE_SENS" =
    setequal(unname(covs_analysis), COVS_SCORE_SENS),
  "truth_sensitivity must have one row per cell" =
    setequal(truth_sensitivity$scenario, superpop_results$scenario) &&
    !anyDuplicated(truth_sensitivity$scenario),
  "exactly one control cell" = sum(superpop_results$control) == 1
)


# ── Cells ─────────────────────────────────────────────────────────────────────
cells <- superpop_results |> select(scenario, beta_B, rho_B_TS, control)

# Each cell carries its own registry and trial pool: workers receive only the
# pools of the cell they run.
cell_configs <- lapply(seq_len(nrow(superpop_results)), \(i) list(
  beta_B        = superpop_results$beta_B[i],
  rho_B_TS      = superpop_results$rho_B_TS[i],
  control       = superpop_results$control[i],
  registry_pool = superpop_results$registry_pool[[i]],
  trial_pool    = superpop_results$trial_pool[[i]]
)) |>
  setNames(superpop_results$scenario)

rm(superpop_results)


# ── Checkpoints ───────────────────────────────────────────────────────────────
run_params <- list(N_BOOT = N_BOOT, N_TRIAL = N_TRIAL, TAU = TAU,
                   RATIO_CONTROLS = RATIO_CONTROLS, N_REGISTRY = N_REGISTRY)

ckpt_dir <- here::here("results", sprintf("checkpoints_sensitivity_%d", N_BOOT))
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

ckpt_file <- function(sc_name) file.path(ckpt_dir, paste0(sc_name, ".rds"))

run_one_cell <- function(sc_name, cfg) {
  file        <- ckpt_file(sc_name)
  truth_cell  <- truth_sensitivity |> dplyr::filter(scenario == sc_name)
  fingerprint <- list(params = run_params, beta_B = cfg$beta_B,
                      rho_B_TS = cfg$rho_B_TS, truth = truth_cell)
  
  if (file.exists(file)) {
    ck <- readRDS(file)
    if (is.list(ck) && identical(ck$fingerprint, fingerprint)) {
      message(sprintf("[cached] %s", sc_name))
      return(data.frame(scenario = sc_name, status = "cached", error_msg = NA_character_))
    }
    message(sprintf("[rerun] %s: checkpoint built with different parameters", sc_name))
  }
  
  # run_scenario_analysis() reads `trial_pool_scenarios` and `seeds` from the
  # global environment, but future's automatic export only sees the objects
  # named in this function's code: both are placed explicitly in the global
  # environment of the worker.
  assign("trial_pool_scenarios", setNames(list(cfg$trial_pool), sc_name), envir = globalenv())
  assign("seeds", seeds, envir = globalenv())
  
  tryCatch({
    res <- run_scenario_analysis(
      label                         = sc_name,
      scenario_name                 = sc_name,
      trial_pool_name               = sc_name,
      sc                            = scenario1,
      B_setting                     = "all",
      include_B_in_score            = FALSE,
      include_unneccessary_in_score = FALSE,
      registry_pool                 = cfg$registry_pool,
      truth_scenarios               = truth_cell,
      n_reg_values                  = N_REGISTRY,
      n_boot                        = N_BOOT,
      n_trial                       = N_TRIAL,
      ratio_controls                = RATIO_CONTROLS,
      tau                           = TAU
    ) |>
      dplyr::mutate(beta_B = cfg$beta_B, rho_B_TS = cfg$rho_B_TS, control = cfg$control)
    
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
# Each cell carries two pools (registry and trial): raise the export limit
options(future.globals.maxSize = 4 * 1024^3)

n_workers <- min(length(cell_configs), max(1L, parallelly::availableCores() - 1L))
plan(multisession, workers = n_workers)
message(sprintf("Running %d cells on %d workers.", length(cell_configs), n_workers))

run_log <- furrr::future_imap(
  cell_configs,
  \(cfg, sc_name) run_one_cell(sc_name, cfg),
  .options = furrr::furrr_options(seed = NULL)
) |> dplyr::bind_rows()

plan(sequential)
print(run_log)

if (any(run_log$status == "error")) {
  stop("Cell(s) failed: ",
       paste(run_log$scenario[run_log$status == "error"], collapse = ", "),
       ". Nothing saved.")
}


# ── Collect results from the checkpoints of the current cells ─────────────────
res_sens <- purrr::map(names(cell_configs), \(s) readRDS(ckpt_file(s))$res) |>
  dplyr::bind_rows()

# Crashed iterations (code failures) are kept apart; match_failed (degenerate
# matching) stays in the results as a reportable outcome.
res_sens_errors <- res_sens |> dplyr::filter(status == "error")
res_sens        <- res_sens |> dplyr::filter(status == "ok")

if (nrow(res_sens_errors) > 0) {
  warning(sprintf("%d iteration(s) crashed; see res_sens_errors.", nrow(res_sens_errors)))
  print(dplyr::count(res_sens_errors, scenario, registry_n, error_msg))
}

stopifnot(
  "effect_setting missing in some rows" =
    all(res_sens$effect_setting %in% c("null", "alternative"))
)

iter_check <- res_sens |>
  dplyr::distinct(scenario, registry_n, effect_setting, iter) |>
  dplyr::count(scenario, registry_n, effect_setting, name = "n_iter") |>
  dplyr::filter(n_iter != N_BOOT)
if (nrow(iter_check) > 0) {
  warning("Some cells have fewer than N_BOOT iterations:")
  print(iter_check)
}

res_sens_alt  <- res_sens |> dplyr::filter(effect_setting == "alternative")
res_sens_null <- res_sens |> dplyr::filter(effect_setting == "null")


# ── Save ──────────────────────────────────────────────────────────────────────
output_file <- here::here(sprintf("results/results_sensitivity_raw_%d.RData", N_BOOT))

save(
  res_sens_alt, res_sens_null, res_sens_errors,
  truth_sensitivity, cells,
  N_BOOT, N_TRIAL, TAU, RATIO_CONTROLS, N_REGISTRY, ALPHA,
  file = output_file
)

message(sprintf("Done. Results saved to %s", output_file))