#senza output
run_scenario_analysis <- function(label,
                                  scenario_name,
                                  trial_pool_name,
                                  sc,
                                  B_setting,
                                  include_B_in_score,
                                  include_unneccessary_in_score = FALSE,
                                  registry_pool,
                                  truth_scenarios,
                                  n_reg_values = c(75,100, 125, 175, 250, 500, 1000),
                                  n_boot = N_BOOT,
                                  n_trial = N_TRIAL,
                                  ratio_controls = RATIO_CONTROLS,
                                  tau = TAU,
                                  effect_setting = c("alternative", "null")) {
  
  effect_setting <- match.arg(effect_setting)
  
  truth <- truth_scenarios |>
    dplyr::filter(scenario == scenario_name)
  
  if (nrow(truth) != 1) {
    stop("truth_scenarios must contain exactly one row for scenario: ", scenario_name)
  }
  
  covs <- make_covariates(
    sc = sc,
    include_B = include_B_in_score,
    include_unneccessary = include_unneccessary_in_score
  )
  
  trial_pool <- trial_pool_scenarios[[trial_pool_name]]
  
  trial_params <- sc |> dplyr::filter(parameters == "trial")
  trial_params$n_trial <- n_trial
  
  res_all <- list()
  
  for (iter in seq_len(n_boot)) {
    
    if (iter %% 500 == 0) {
      message(sprintf("[%s] %s | %s | iter %d / %d",
                      format(Sys.time(), "%H:%M:%S"),
                      scenario_name, effect_setting, iter, n_boot))
    }
    
    nested_registries <- make_nested_registries(
      registry_pool = registry_pool,
      n_reg_values  = n_reg_values,
      seed          = seeds[iter]
    )
    
    for (namereg in names(nested_registries)) {
      
      registry <- nested_registries[[namereg]]
      
      models <- train_models(
        registry   = registry,
        covariates = unname(covs)
      )
      
      res_iter <- tryCatch(one_resample2(
        iter           = iter,
        prog_model     = models,
        seeds          = seeds,
        trial_pool     = trial_pool,
        registry       = registry,
        B              = B_setting,
        n_trial        = trial_params$n_trial,
        truth          = truth,
        ratio_controls = ratio_controls,
        covs           = covs,
        tau            = tau,
        effect_setting = effect_setting
      ) |> dplyr::mutate(status = "ok"),
      error = function(e) {
        data.frame(
          iter      = iter,
          seed      = seeds[iter],
          status    = "error",
          error_msg = conditionMessage(e)
        )
      }
      )
      
      res_iter <- res_iter |>
        dplyr::mutate(
          scenario           = scenario_name,
          analysis_label     = label,
          registry_n         = namereg,
          B_setting          = B_setting,
          include_B_in_score = include_B_in_score
          #,effect_setting     = effect_setting
        )
      
      res_all[[paste(iter, namereg, sep = "_")]] <- res_iter
    }
  }
  
  dplyr::bind_rows(res_all)  # restituisce solo res_boot
}