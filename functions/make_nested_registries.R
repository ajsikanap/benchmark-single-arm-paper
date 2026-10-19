make_nested_registries <- function(registry_pool,
                                   n_reg_values = c(75,100, 125, 175, 250, 500, 1000),
                                   seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  max_n <- max(n_reg_values)
  
  if (nrow(registry_pool) < max_n) {
    stop("registry_pool has fewer rows than max(n_reg_values).")
  }
  
  registry_big <- registry_pool |>
    dplyr::slice_sample(n = max_n, replace = FALSE)
  
  out <- purrr::map(
    n_reg_values,
    ~ registry_big |> dplyr::slice_head(n = .x)
  )
  
  names(out) <- paste0(n_reg_values)
  
  out
}