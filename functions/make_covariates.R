make_covariates <- function(sc, include_B = TRUE, include_unneccessary=FALSE) {
  
  covs <- sc |>
    dplyr::filter(parameters == "trial") |>
    dplyr::select(
      M, age, sex, resection, trt_type, WHOm2,
      # localization_0,
      # localization_1,
      # localization_2,
      # localization_3,
      # localization_4,
      # localization_5,
      multifocal,
      pTERT, tumor_size, EGFR,
      resection,
      g7l10, 
      unneccessary
    )
  
  covs_vector <- names(covs)[which(covs[1, ] == 1)]
  
  if (include_B) {
    covs_vector <- c("B", covs_vector)
  }
  if (include_unneccessary) {
    covs_vector <- c("unneccessary", covs_vector)
  }
  
  # one_resample2/run_matching expects named covariates
  stats::setNames(covs_vector, covs_vector)
}

