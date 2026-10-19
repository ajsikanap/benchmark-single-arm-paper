
# this is not needed anymore. Error due to summarize and then mutate, but summarize does
# contain only n_iter.
summarise_bias_decomposition <- function(res_boot,
                                         group_vars = c("method"),
                                         return_flextable = FALSE) {
  mcse_bias <- function(b) {
    n <- sum(!is.na(b)); if (n < 2) return(NA_real_)
    sd(b, na.rm = TRUE) / sqrt(n)
  }
  mcse_rmse <- function(b) {
    n  <- sum(!is.na(b)); rm <- sqrt(mean(b^2, na.rm = TRUE))
    if (is.na(rm) || rm == 0 || n < 2) return(NA_real_)
    sd(b^2, na.rm = TRUE) / (2 * rm * sqrt(n))
  }
  mcse_prop <- function(ind) {
    ind <- ind[!is.na(ind)]; n <- length(ind); if (n == 0) return(NA_real_)
    p <- mean(ind); sqrt(p * (1 - p) / n)
  }
  
  out <- res_boot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      
      n_iter = dplyr::n(),
      
      
      .groups = "drop"
    ) |>
    # Check algebrico: deve essere ~0 in ogni scenario
    # Se diverso da 0 indica un bug nella decomposizione
    dplyr::mutate(
      check_deltaS_matched  = bias_deltaS_matched   - reconstructed_bias_deltaS_matched,
      check_deltaS_naive    = bias_deltaS_naive     - reconstructed_bias_deltaS_naive,
      check_deltaRMST_matched = bias_deltaRMST_matched - reconstructed_bias_deltaRMST_matched,
      check_deltaRMST_naive   = bias_deltaRMST_naive   - reconstructed_bias_deltaRMST_naive
    )
  
  if (return_flextable) flextable::flextable(out) else out
}

