summarise_endpoint_components <- function(res_boot,
                                          group_vars = c("method"),
                                          return_flextable = FALSE) {
  
  out <- res_boot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      
      n_iter = dplyr::n(),
      
      
      
      .groups = "drop"
    )
  
  if (return_flextable) flextable::flextable(out) else out
}


summarise_endpoint_components2 <- function(res_boot,
                                           group_vars = c("method"),
                                           formatted = TRUE,
                                           return_flextable = FALSE) {
  
  endpoint_vars <- c(
    # Survival latent sampled values
    "trial_true_z1",
    "trial_true_z0",
    "registry_true",
    "matched_true_z0",
    
    # RMST latent sampled values
    "rmst_trial_true_z1",
    "rmst_trial_true_z0",
    "rmst_registry_true",
    "rmst_matched_true_z0",
    
    # Survival observed/KM estimates
    "trial_z1",
    "trial_z0",
    "matched",
    "registry",
    
    # RMST observed/KM estimates
    "rmst_trial_z1",
    "rmst_trial_z0",
    "rmst_matched",
    "rmst_registry"
  )
  
  # Input checks
  missing_group_vars <- setdiff(group_vars, names(res_boot))
  if (length(missing_group_vars) > 0) {
    stop("These group_vars are missing from res_boot: ",
         paste(missing_group_vars, collapse = ", "))
  }
  
  missing_endpoint_vars <- setdiff(endpoint_vars, names(res_boot))
  if (length(missing_endpoint_vars) > 0) {
    stop("These endpoint variables are missing from res_boot: ",
         paste(missing_endpoint_vars, collapse = ", "))
  }
  
  format_mean_iqr_sd <- function(x, digits = 3) {
    m  <- mean(x,                         na.rm = TRUE)
    s  <- sd(x,                           na.rm = TRUE)
    q1 <- as.numeric(quantile(x, 0.25,   na.rm = TRUE))
    q3 <- as.numeric(quantile(x, 0.75,   na.rm = TRUE))
    sprintf(
      paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f] (%.", digits, "f)"),
      m, q1, q3, s
    )
  }
  
  # Benchmark alignment differences
  # Convenzione: benchmark - trial_cf
  # positivo = benchmark piu favorevole del contraffattuale (benchmark ottimistico)
  # negativo = benchmark meno favorevole del contraffattuale (benchmark pessimistico, caso S4)
  out <- res_boot |>
    dplyr::mutate(
      
      # Survival latent benchmark alignment
      diff_naive_bench_true   = registry_true    - trial_true_z0,
      diff_matched_bench_true = matched_true_z0  - trial_true_z0,
      
      # RMST latent benchmark alignment
      rmst_diff_naive_bench_true   = rmst_registry_true   - rmst_trial_true_z0,
      rmst_diff_matched_bench_true = rmst_matched_true_z0 - rmst_trial_true_z0,
      
      # Survival estimated benchmark alignment
      # Nota: trial_z0 e' stimato sui tempi latenti — quantita' diagnostica,
      # non osservabile in un trial reale
      diff_naive_bench_est   = registry - trial_z0,
      diff_matched_bench_est = matched  - trial_z0,
      
      # RMST estimated benchmark alignment
      rmst_diff_naive_bench_est   = rmst_registry - rmst_trial_z0,
      rmst_diff_matched_bench_est = rmst_matched  - rmst_trial_z0
    )
  
  summary_vars <- c(
    endpoint_vars,
    "diff_naive_bench_true",
    "diff_matched_bench_true",
    "rmst_diff_naive_bench_true",
    "rmst_diff_matched_bench_true",
    "diff_naive_bench_est",
    "diff_matched_bench_est",
    "rmst_diff_naive_bench_est",
    "rmst_diff_matched_bench_est"
  )
  
  if (formatted) {
    out <- out |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        n_iter = dplyr::n(),
        dplyr::across(
          dplyr::all_of(summary_vars),
          ~ format_mean_iqr_sd(.x),
          .names = "{.col}"
        ),
        .groups = "drop"
      )
  } else {
    out <- out |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        n_iter = dplyr::n(),
        dplyr::across(
          dplyr::all_of(summary_vars),
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd   = ~ sd(.x,   na.rm = TRUE),
            q25  = ~ as.numeric(quantile(.x, 0.25, na.rm = TRUE)),
            q75  = ~ as.numeric(quantile(.x, 0.75, na.rm = TRUE))
          ),
          .names = "{.col}_{.fn}"
        ),
        .groups = "drop"
      )
  }
  
  if (return_flextable) {
    flextable::flextable(out)
  } else {
    out
  }
}

summarise_endpoint_components3 <- function(res_boot,
                                           group_vars = c("method"),
                                           return_flextable = FALSE) {
  
  endpoint_vars <- c(
    ## Survival latent sampled values
    "trial_true_z1",
    "trial_true_z0",
    "registry_true",
    "matched_true_z0",
    
    ## RMST latent sampled values
    "rmst_trial_true_z1",
    "rmst_trial_true_z0",
    "rmst_registry_true",
    "rmst_matched_true_z0",
    
    ## Survival observed/KM estimates
    "trial_z1",
    "trial_z0",
    "matched",
    "registry",
    
    ## RMST observed/KM estimates
    "rmst_trial_z1",
    "rmst_trial_z0",
    "rmst_matched",
    "rmst_registry"
  )
  
  missing_group_vars <- setdiff(group_vars, names(res_boot))
  if (length(missing_group_vars) > 0) {
    stop(
      "These group_vars are missing from res_boot: ",
      paste(missing_group_vars, collapse = ", ")
    )
  }
  
  missing_endpoint_vars <- setdiff(endpoint_vars, names(res_boot))
  if (length(missing_endpoint_vars) > 0) {
    stop(
      "These endpoint variables are missing from res_boot: ",
      paste(missing_endpoint_vars, collapse = ", ")
    )
  }
  
  format_mean_minmax_sd <- function(x, digits = 3) {
    m  <- mean(x, na.rm = TRUE)
    s  <- stats::sd(x, na.rm = TRUE)
    q1 <- quantile(x, probs=0.25, na.rm = TRUE)
    q3 <- quantile(x, probs=0.75, na.rm = TRUE)
    
    sprintf(
      paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f] (%.", digits, "f)"),
      m, q1, q3, s
    )
  }
  
  out <- res_boot |>
    dplyr::mutate(
      ## Survival latent benchmark differences
      diff_naive_bench_true =
        trial_true_z0 - registry_true,
      
      diff_matched_bench_true =
        trial_true_z0 - matched_true_z0,
      
      ## RMST latent benchmark differences
      rmst_diff_naive_bench_true =
        rmst_trial_true_z0 - rmst_registry_true,
      
      rmst_diff_matched_bench_true =
        rmst_trial_true_z0 - rmst_matched_true_z0,
      
      ## Survival observed/KM benchmark differences
      diff_naive_bench_est =
        trial_z0 - registry,
      
      diff_matched_bench_est =
        trial_z0 - matched,
      
      ## RMST observed/KM benchmark differences
      rmst_diff_naive_bench_est =
        rmst_trial_z0 - rmst_registry,
      
      rmst_diff_matched_bench_est =
        rmst_trial_z0 - rmst_matched
    )
  
  summary_vars <- c(
    endpoint_vars,
    "diff_naive_bench_true",
    "diff_matched_bench_true",
    "rmst_diff_naive_bench_true",
    "rmst_diff_matched_bench_true",
    "diff_naive_bench_est",
    "diff_matched_bench_est",
    "rmst_diff_naive_bench_est",
    "rmst_diff_matched_bench_est"
  )
  
  out <- out |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_vars))
    ) |>
    dplyr::summarise(
      n_iter = dplyr::n(),
      dplyr::across(
        dplyr::all_of(summary_vars),
        ~ format_mean_minmax_sd(.x),
        .names = "{.col}"
      ),
      .groups = "drop"
    )
  
  if (return_flextable) {
    flextable::flextable(out)
  } else {
    out
  }
}
