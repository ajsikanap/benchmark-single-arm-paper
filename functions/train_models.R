train_models_non_penalized <- function(registry, covariates) {
  
  ## prognostic model
  prog_formula <- as.formula(
    paste("Surv(time, event) ~", paste(covariates, collapse = " + "))
  )
  prog_model <- coxph(prog_formula, data = registry)
  
  ## propensity model (registry vs trial pool logic is WRONG here)
  ## this must be trained later, on (trial + registry) by design
  list(prog_model = prog_model)
}

# install.packages("coxphf")
train_models <- function(registry, covariates) {
  
  prog_formula_standard <- as.formula(
    paste("Surv(time, event) ~", paste(covariates, collapse = " + "))
  )
  
  # ridge() va nella formula, non come argomento di coxph
  prog_formula_ridge <- as.formula(
    paste("Surv(time, event) ~ ridge(",
          paste(covariates, collapse = ", "),
          ", theta = 0.1, scale = TRUE)")
  )
  
  prog_model <- tryCatch(
    #survival::coxph(prog_formula_ridge, data = registry),
    survival::coxph(prog_formula_standard, data = registry),
    error = function(e) {
      warning("Ridge Cox failed, falling back to standard coxph: ",
              conditionMessage(e))
      survival::coxph(prog_formula_standard, data = registry)
    }
  )
  
  list(prog_model = prog_model)
}
