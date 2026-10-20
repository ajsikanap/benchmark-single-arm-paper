# ============================================================================
# 08_run_sensitivity.R
#
# Run della sensitivity analysis S5 sulla griglia (beta_B, rho_B_TS).
# Riusa run_scenario_analysis() (functions/run_scenario_analysis0.R) e
# one_resample2() = one_resample2v2 (functions/one_resample2v2.R), gia'
# validate per S0-S4. Nessuna nuova funzione di analisi qui: questo script
# produce SOLO i res_boot grezzi (un data.frame per cella, con la colonna
# effect_setting in {"null","alternative"}).
#

# QUi chiamo una sola volta per cella run_analysis; 
# il filtro per effect_setting va fatto nello
# script di analisi (08_analyse_sensitivity.Rmd), separatamente per T1E
# (effect_setting=="null") e power/bias (effect_setting=="alternative").
# ============================================================================


# ── Functions and inputs ──────────────────────────────────────────────────────
source(here::here("functions/load_functions.R"))
source(here::here("config/parameters.R"))

N_REGISTRY <- sort(unique(c(N_REGISTRY_SENS, max(N_REG_VALUES))))

scenario1 <- readxl::read_excel(here::here(PATH_DGP_XLSX), sheet = "S1")

seeds <- read.table(here::here(SEED_FILE), header = FALSE)$V1

load(here::here("data/superpop_results_sensitivity.rda"))   # superpop_results
load(here::here("data/truth_sensitivity.rda"))              # truth_sensitivity



# # ── Parallelizzazione ─────────────────────────────────────────────────────────
# Usa tutti i core disponibili meno uno (lascia un core libero per il sistema).
# Su un server con molti core puoi alzare il numero o togliere il -1.
n_cores <- max(1L, parallelly::availableCores() - 1L)
plan(multisession, workers = n_cores)
message(sprintf("Activate parallel computing: %d worker.", n_cores))

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
N_BOOT         <- 10
N_REGISTRY     <- 250  # singola dimensione di registry per S5

# ---------------------------------------------------------------------------
# Funzioni e dati di base (stesse di S0-S4)
# ---------------------------------------------------------------------------

s1 <- scenario1
reg_params <- s1 |> dplyr::filter(parameters == "registry")

# ---------------------------------------------------------------------------
# Pool registry/trial per cella, con colonna S richiesta da score_data()
# ---------------------------------------------------------------------------
sensitivity_pools <- superpop_results |>
  dplyr::mutate(
    registry_pool = purrr::map(superpop_obj, "registry_pool"),
    trial_pool    = purrr::map(superpop_obj, "trial_pool"),
    checks        = purrr::map(superpop_obj, "checks"),
    scenario_id   = paste0(
      "betaB_", beta_B,
      "_rho_", rho_B_TS
    )
  ) |>
  dplyr::select(
    scenario_id,
    beta_B,
    rho_B_TS,
    registry_pool,
    trial_pool,
    checks
  )

# ---------------------------------------------------------------------------
# Griglia S5: 3 (beta_B) x 5 (rho_B_TS) = 15 celle
# ---------------------------------------------------------------------------
main_grid <- tidyr::expand_grid(
  beta_B   = c(-1, 0, 1),
  rho_B_TS = c(-0.75, -0.5, 0, 0.5, 0.75)
)

sensitivity_pools <- sensitivity_pools |>
  dplyr::semi_join(
    main_grid,
    by = c("beta_B", "rho_B_TS")
  )

stopifnot(nrow(sensitivity_pools) == nrow(main_grid))

# ---------------------------------------------------------------------------
# Controllo di coerenza: le covariate usate per la "verita'" in
# truth_sensitivity (covs_score_sensitivity, senza B) devono coincidere con
# quelle che run_scenario_analysis() costruira' internamente via
# make_covariates(sc = s1, include_B = FALSE, include_unneccessary = FALSE).
# Se questo controllo fallisce, beta_adj_prog_star/beta_adj_prop_star in
# truth_sensitivity NON sono le "verita'" corrette per il modello di matching
# effettivamente usato qui.
# ---------------------------------------------------------------------------
covs_check <- make_covariates(
  sc = s1,
  include_B = FALSE,
  include_unneccessary = FALSE
)

covs_score_sensitivity <- c(
  "M", "age", "sex", "WHOm2",
  "resection", "trt_type",
  # "localization_0",
  # "localization_1",
  # "localization_2",
  # "localization_3",
  # "localization_4",
  # "localization_5",
  "multifocal",
  "pTERT", "tumor_size",
  "EGFR", "g7l10"
)

if (!setequal(unname(covs_check), covs_score_sensitivity)) {
  stop(
    "Le covariate di make_covariates(s1, include_B=FALSE) non coincidono ",
    "con covs_score_sensitivity usate per calcolare truth_sensitivity.\n",
    "make_covariates(): ", paste(unname(covs_check), collapse = ", "), "\n",
    "covs_score_sensitivity: ", paste(covs_score_sensitivity, collapse = ", ")
  )
}

# ---------------------------------------------------------------------------
# Una sola run per cella (NIENTE doppia chiamata alt/null)
# ---------------------------------------------------------------------------
run_sensitivity_cell <- function(row_id) {
  
  cell <- sensitivity_pools |>
    dplyr::slice(row_id)
  
  scenario_id <- cell$scenario_id
  
  registry_pool_cell <- cell$registry_pool[[1]] |> dplyr::mutate(S = 0)
  trial_pool_cell    <- cell$trial_pool[[1]]    |> dplyr::mutate(S = 1)
  
  truth_cell <- truth_sensitivity |>
    dplyr::filter(scenario == scenario_id)
  
  # run_scenario_analysis() si aspetta trial_pool_scenarios in .GlobalEnv
  assign(
    "trial_pool_scenarios",
    list(sensitivity = trial_pool_cell),
    envir = .GlobalEnv
  )
  
  # effect_setting e' richiesto dalla firma (match.arg) ma viene ignorato da
  # one_resample2v2: il res_boot restituito contiene GIA' sia le righe
  # "null" sia le righe "alternative".
  run_scenario_analysis(
    label                          = scenario_id,
    scenario_name                  = scenario_id,
    trial_pool_name                = "sensitivity",
    sc                             = s1,
    B_setting                      = "all",
    include_B_in_score             = FALSE,
    include_unneccessary_in_score  = FALSE,
    registry_pool                  = registry_pool_cell,
    truth_scenarios                = truth_cell,
    n_reg_values                   = N_REGISTRY,
    n_boot                         = N_BOOT,
    n_trial                        = N_TRIAL,
    ratio_controls                 = RATIO_CONTROLS,
    tau                            = TAU,
    effect_setting                 = "alternative"
  )
}

# ---------------------------------------------------------------------------
# Esecuzione su tutte le celle
# ---------------------------------------------------------------------------
n_cells <- nrow(sensitivity_pools)

results_list_sensitivity <- vector("list", n_cells)
names(results_list_sensitivity) <- sensitivity_pools$scenario_id

for (i in seq_len(n_cells)) {
  
  # message outputted for kkeping track of sim
  message(sprintf(
    "[%s] cella %d/%d: %s (beta_B=%.2f, rho_B_TS=%.2f)",
    format(Sys.time(), "%H:%M:%S"),
    i, n_cells,
    sensitivity_pools$scenario_id[i],
    sensitivity_pools$beta_B[i],
    sensitivity_pools$rho_B_TS[i]
  ))
  
  results_list_sensitivity[[i]] <- run_sensitivity_cell(row_id = i)
}

# ---------------------------------------------------------------------------
# Salvataggio
#
# Ogni elemento di results_list_sensitivity e' un res_boot piatto
# (n_boot x 2 metodi x 2 effect_setting righe), con colonne:
# scenario, registry_n="250", method in {prog,prop},
# effect_setting in {null, alternative}, e tutte le colonne di output di
# one_resample2v2 (bias_*, err_*, *_p_one_sided_*, match_diag, ecc.).
# ---------------------------------------------------------------------------
output_file <-here::here(
  sprintf("results/results_sensitivity_betaB_rho_TS_registry250_%d.Rdata", N_BOOT)
)


save(
  results_list_sensitivity,
  sensitivity_pools,
  truth_sensitivity,
  N_BOOT, N_TRIAL, TAU, RATIO_CONTROLS, N_REGISTRY,
  file =output_file
)



message(sprintf("Results saved to %s",output_file))
