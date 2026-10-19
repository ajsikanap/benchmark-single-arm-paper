one_sample_logrank_km <- function(time, event, hist_time, hist_event, alpha=0.1) {
  
  fit_hist <- survival::survfit(
    survival::Surv(hist_time, hist_event) ~ 1
  )
  
  hist_df <- data.frame(
    time = fit_hist$time,
    surv = fit_hist$surv
  )
  
  hist_df <- hist_df |>
    dplyr::filter(surv > 0)
  
  hist_df$H0 <- -log(hist_df$surv)
  
  H0_fun <- stats::stepfun(
    x = hist_df$time,
    y = c(0, hist_df$H0)
  )
  
  O <- sum(event)
  E <- sum(H0_fun(time))
  V <- E
  
  Z <- (O - E) / sqrt(V)
  
  p_two_sided <- 2 * (1 - pnorm(abs(Z)))
  
  # If treatment is expected to improve survival, fewer events than expected:
  # O < E, hence Z < 0.
  p_one_sided_better <- pnorm(Z)
  
  data.frame(
    one_sample_lr_O = O,
    one_sample_lr_E = E,
    one_sample_lr_V = V,
    one_sample_lr_Z = Z,
    one_sample_lr_p_two_sided = p_two_sided,
    one_sample_lr_p_one_sided_better = p_one_sided_better,
    one_sample_lr_reject = as.integer(p_one_sided_better < alpha)
  )
}