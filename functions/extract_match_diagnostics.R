

extract_match_diagnostics <- function(m.out, smd_dt = NULL) {
  
  sm <- summary(m.out)
  
  score_before <- sm$sum.all
  score_after  <- sm$sum.matched
  
  out <- tibble::tibble(
    score_smd_before = abs(score_before["distance", "Std. Mean Diff."]),
    score_smd_matched = abs(score_after["distance", "Std. Mean Diff."]),
    
    score_var_ratio_before = score_before["distance", "Var. Ratio"],
    score_var_ratio_matched = score_after["distance", "Var. Ratio"],
    
    score_pair_distance = score_after["distance", "Std. Pair Dist."],
    
    n_control_all = sm$nn["All", "Control"],
    n_trial_all = sm$nn["All", "Treated"],
    n_control_matched = sm$nn["Matched", "Control"],
    # n_trial_matched = sm$nn["Matched", "Treated"],
    n_control_ess = sm$nn["Matched (ESS)", "Control"],
    n_control_discarded = sm$nn["Discarded", "Control"],
    # n_trial_discarded = sm$nn["Discarded", "Treated"],
    
    prop_control_discarded = n_control_discarded / n_control_all
    # prop_trial_matched = n_trial_matched / n_trial_all
  ) |>
    dplyr::mutate(
      score_smd_reduction_abs = score_smd_before - score_smd_matched,
      score_smd_reduction_rel = score_smd_reduction_abs / score_smd_before
    )
  
  if (!is.null(smd_dt)) {
    
    cov_smd <- smd_dt$Balance |>
      tibble::rownames_to_column("covariate") |>
      dplyr::filter(!covariate %in% c("distance", "score")) |>
      dplyr::mutate(abs_smd = abs(Diff.Adj))
    
    out <- out |>
      dplyr::mutate(
        smd_cov_mean_matched = mean(cov_smd$abs_smd, na.rm = TRUE),
        smd_cov_max_matched = max(cov_smd$abs_smd, na.rm = TRUE),
        smd_cov_min_matched = min(cov_smd$abs_smd, na.rm = TRUE),
        #proportion of cov with smd over 0.1
        smd_cov_prop_gt_0.1 = mean(cov_smd$abs_smd > 0.1, na.rm = TRUE),
        #proportion of cov with smd over 0.2
        smd_cov_prop_gt_0.2 = mean(cov_smd$abs_smd > 0.2, na.rm = TRUE)
      )
  }
  
  out
}
