get_surv_12 <- function(fit, time = 12) {
  if (is.null(fit)) return(NA_real_)
  out <- tryCatch(
    summary(fit, times = time)$surv,
    error = function(e) NA_real_
  )
  if (length(out) == 0) NA_real_ else out
}

# 4. Endpoint assessment
get_rmst_from_km <- function(fit, tau = 12) {
  if (is.null(fit)) return(NA_real_)
  ss <- summary(fit)
  keep <- ss$time < tau
  times <- c(0, ss$time[keep], tau)
  survs <- c(1, ss$surv[keep])
  sum(diff(times) * survs)
}