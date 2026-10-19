

estimate_treatment <- function(matched_data, beta_true, method = c("source_only", "source_prog_score", "source_score")) {
  
  method <- match.arg(method)
  
  # Ensure strata_id exists if needed
  if(method == "source_strata" & !"strata_id" %in% colnames(matched_data)) {
    matched_data$strata_id <- matched_data$subclass
  }
  
  # Define formula based on method
  formula <- switch(
    method,
    source_only   = Surv(time, event) ~ source,
    source_prog_score = Surv(time, event) ~ source + prog_score,
    source_score  = Surv(time, event) ~ source + score
  )
  
  na_row <- data.frame(beta_hat = NA_real_, se_hat = NA_real_, ci_low = NA_real_,
                       ci_high = NA_real_, cover = NA, ic_len = NA_real_,
                       p_val = NA_real_, sign_pval = NA)
  # Fit Cox model
  fit <- tryCatch(coxph(formula, data = matched_data), error = function(e) NULL)
  if (is.null(fit)) return(na_row)
  cf <- coef(fit)
  if (!"source" %in% names(cf) || is.na(cf[["source"]])) return(na_row)
  
  beta_hat <- unname(cf[["source"]])
  vc <- vcov(fit)
  se_hat  <- if ("source" %in% rownames(vc)) sqrt(vc["source","source"]) else NA_real_
  # 
  fit <- coxph(formula, data = matched_data)
  
  beta_hat <- coef(fit)["source"]
  se_hat   <- sqrt(vcov(fit)["source","source"])
  
  ci_low  <- beta_hat - qnorm(0.975) * se_hat
  ci_high <- beta_hat + qnorm(0.975) * se_hat
  
  sc <- summary(fit)$coefficients
  p_val <- if ("source" %in% rownames(sc)) sc["source","Pr(>|z|)"] else NA_real_

  
  data.frame(beta_hat, se_hat, ci_low, ci_high,
             cover = beta_true >= ci_low & beta_true <= ci_high,
             ic_len = ci_high - ci_low, 
             p_val = p_val,
             sign_pval = (beta_hat / se_hat) < qnorm(0.10))
  
}