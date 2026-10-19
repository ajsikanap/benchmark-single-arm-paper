summarise_balance_results <- function(res_boot,
                                      group_vars = c("registry_n", "method"),
                                      digits = 2,
                                      return_flextable = FALSE) {
  
  fmt_mean_iqr_sd <- function(x, digits = 2) {
    m   <- mean(x, na.rm = TRUE)
    q25 <- stats::quantile(x, 0.25, na.rm = TRUE)
    q75 <- stats::quantile(x, 0.75, na.rm = TRUE)
    s   <- stats::sd(x, na.rm = TRUE)
    
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f); SD=%.", digits, "f"),
      m, q25, q75, s
    )
  }
  
  metric_vars <- c(
    "score_smd_before",
    "score_smd_matched",
    "score_pair_distance",
    "smd_cov_mean_matched",
    "smd_cov_prop_gt_0.1",
    "smd_cov_prop_gt_0.2",
    "n_control_ess",
    "n_control_discarded",
    "prop_control_discarded"
  )
  
  out <- res_boot |>
    dplyr::select(
      dplyr::all_of(group_vars),
      dplyr::all_of(metric_vars)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      
      `SMD score unmatched`       = fmt_mean_iqr_sd(score_smd_before,            digits),
      `SMD score matched`         = fmt_mean_iqr_sd(score_smd_matched,           digits),
      `Mean pair distance`        = fmt_mean_iqr_sd(score_pair_distance,         digits),
      `Mean covariate SMD matched`= fmt_mean_iqr_sd(smd_cov_mean_matched,        digits),
      `% covariates SMD > 0.10`   = fmt_mean_iqr_sd(100 * smd_cov_prop_gt_0.1,  digits),
      `% covariates SMD > 0.20`   = fmt_mean_iqr_sd(100 * smd_cov_prop_gt_0.2,  digits),
      `ESS matched controls`      = fmt_mean_iqr_sd(n_control_ess,               digits),
      `Controls discarded (n)`    = fmt_mean_iqr_sd(n_control_discarded,         digits),
      `Controls discarded (%)`    = fmt_mean_iqr_sd(100 * prop_control_discarded,digits),
      
      .groups = "drop"
    ) |>
    dplyr::arrange(registry_n, method)
  
  if (return_flextable) {
    out <- out |>
      flextable::flextable() |>
      flextable::autofit()
  }
  
  out
}