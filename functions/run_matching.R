# score data
score_data <- function(trial, registry, prog_model) {
  
  registry$prog_score <- predict(prog_model[[1]], type = "lp")
  trial$prog_score    <- predict(prog_model[[1]], newdata = trial, type = "lp")
  
  rbind(
    transform(trial |> dplyr::select(-S), source = 1),
    transform(registry  |> dplyr::select(-S), source = 0)
  )
}


run_matching <- function(match_data, covariates,
                         B,
                         method = c("prog", "prop"),
                         ratio_controls = 1) {
  
  method <- match.arg(method)
  
  if(B=="unfiltered"){
    match_data<-match_data # |> filter((B==1 & source==1)|(source==0 & ((B==1)|(B==0)) ))
  }
  
  if(B=="all"){
    match_data<-match_data  |> filter((B==1 & source==1)|(source==0 & ((B==1)|(B==0)) ))
  }
  
  # two specific sensitivity tests
  if(B=="1"){
    match_data<-match_data |> filter((B==1 & source==1)| (source==0 & B==1))
  }
  if(B=="0"){
    match_data<-match_data |> filter((B==1 & source==1)| (source==0 & B==0))
  }
  
  
  if (method == "prop") {
    
    ps_formula <- as.formula(
      paste("source ~", paste(names(covariates), collapse = " + "))
    )
    
    
    # ── Propensity score con Firth per gestire separazione ────────────────────
    ps_model <- tryCatch(
      # logistf::logistf(ps_formula, data = match_data, firth = TRUE,
      #                  control = logistf::logistf.control(maxit = 200)),
      glm(ps_formula, data = match_data, family = binomial()),
      error = function(e) {
        # fallback a glm se logistf fallisce (raro)
        warning("logistf failed, falling back to glm: ", conditionMessage(e))
        glm(ps_formula, data = match_data, family = binomial())
      }
    )
    
    
    match_data$score <- predict(ps_model, type = "response")
    
  } else {
    match_data$score <- match_data$prog_score
  }
  
  m.out <- tryCatch(
    MatchIt::matchit(
    source ~ score,
    data = match_data,
    method = "nearest", replace = FALSE,
    distance = match_data$score,
    ratio = ratio_controls,
    discard = "control"
  ),
  error=function(e) {
    warning(sprintf("matchit failed (iter): %s", conditionMessage(e)))
    NULL
  }
  )
  
  if (is.null(m.out)){
    return(list(matched=NULL, m.out=NULL,
                filtered_matched=match_data, match_failed=TRUE))
  }
  
  list(matched=MatchIt::match.data(m.out),m.out=m.out,filtered_matched=match_data, match_failed=FALSE)
}
