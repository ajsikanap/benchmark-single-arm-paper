make_density_plot <- function(res_boot, title = NULL) {
  
  df_long <- res_boot |>
    dplyr::select(iter, method, trial_z1, trial_z0, matched, registry) |>
    tidyr::pivot_longer(
      cols = c(trial_z1, trial_z0, matched, registry),
      names_to = "source",
      values_to = "value"
    ) |>
    dplyr::mutate(
      source = factor(
        source,
        levels = c("trial_z1", "trial_z0", "matched", "registry"),
        labels = c("Trial Z=1", "Trial Z=0", "Matched control", "Registry")
      ),
      method = factor(method, levels = c("prop", "prog"))
    )
  
  p <- ggplot(df_long, aes(x = value, fill = source)) +
    geom_density(alpha = 0.35) +
    facet_wrap(~ method, ncol = 1) +
    labs(
      x = "Survival at 12 months",
      y = "Density",
      fill = ""
    ) +
    theme_minimal()
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  list(
    df_long = df_long,
    plot = p
  )
}