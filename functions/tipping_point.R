# tipping_point.R
# ---------------------------------------------------------------------------
# delta* (tipping point) per il one-sample log-rank e derivati.
#
# Fondamento (Woolson 1981; Finkelstein et al. 2003):
#   OSLRT = (O - E) / sqrt(E),   E = sum_i Lambda0(X_i)
# Perturbazione proporzionale del benchmark: Lambda0 -> theta * Lambda0,
# quindi E -> theta*E, O invariato (e' del trial). Imponendo Z(theta) = -z_alpha:
#
#   theta* = ( z_alpha + sqrt(z_alpha^2 + 4*O) )^2 / (4*E)
#   delta* = log(theta*)                       [log-hazard ratio del benchmark]
#
# delta* < 0  -> decisione corrente GO; |delta*| = margine
# delta* > 0  -> decisione corrente NO-GO; serve benchmark peggiore di e^{delta*}
#
# Richiede in res_boot le colonne prodotte da one_resample2v2:
#   one_sample_lr_O_matched/_registry, one_sample_lr_E_matched/_registry,
#   e (per la scala sopravvivenza) matched, registry.
# ---------------------------------------------------------------------------


# --- 1. delta* per replicato -----------------------------------------------
add_tipping_point <- function(df, alpha = 0.10) {
  z_a <- stats::qnorm(1 - alpha)
  
  tp <- function(O, E) {
    theta <- (z_a + sqrt(z_a^2 + 4 * O))^2 / (4 * E)
    list(theta = theta, delta = log(theta))
  }
  
  m <- tp(df$one_sample_lr_O_matched,  df$one_sample_lr_E_matched)
  r <- tp(df$one_sample_lr_O_registry, df$one_sample_lr_E_registry)
  
  df$theta_star_matched  <- m$theta
  df$delta_star_matched  <- m$delta
  df$theta_star_registry <- r$theta
  df$delta_star_registry <- r$delta
  df
}


# --- 2. conversione delta* -> scala sopravvivenza + rapporto di fragilita' --
# S_tip(tau) = exp(-theta* * H0(tau)),  H0(tau) = -log(S_bench(tau))
# S_bench(tau) empirico = colonna `matched` (o `registry`), gia' KM a tau.
#
# Rapporto di fragilita' (con segno): quanto l'errore strutturale vero occupa
# della tolleranza concessa da delta*. > 1 => la verita' e' oltre il tipping
# point => il GO e' un artefatto. Calcolabile SOLO in simulazione (usa err_matching_z0).
add_survival_scale <- function(df, alpha = 0.10) {
  if (!"theta_star_matched" %in% names(df)) df <- add_tipping_point(df, alpha)
  
  # benchmark matched
  df$H0_tau_matched <- -log(df$matched)
  df$S_tip_matched  <- exp(-df$theta_star_matched * df$H0_tau_matched)
  df$d_star_matched <- df$S_tip_matched - df$matched          # tolleranza (con segno)
  df$ME_matched     <- df$err_matching_z0                     # errore strutturale vero
  df$fragility_matched <- -df$ME_matched / df$d_star_matched
  
  # benchmark registry (naive)
  df$H0_tau_registry <- -log(df$registry)
  df$S_tip_registry  <- exp(-df$theta_star_registry * df$H0_tau_registry)
  df$d_star_registry <- df$S_tip_registry - df$registry
  # per il registro l'errore strutturale e' selection (trial_true - registry_true),
  # con lo stesso segno di err_matching_z0 (matched_true - trial_true).
  df$ME_registry     <- -df$err_registry_selection_z0
  df$fragility_registry <- -df$ME_registry / df$d_star_registry
  
  df
}


# --- 3. aggregato per gruppo (stile summarise_matching_error) ---------------
# Mediana + IQR di delta* e della scala sopravvivenza, quota di replicati
# fragili (>1) e frazione go, per scenario x registry_n x method.
summarise_tipping_point <- function(res_boot,
                                    group_vars = c("scenario", "registry_n", "method"),
                                    alpha = 0.10,
                                    return_flextable = FALSE) {
  
  if (!"delta_star_matched" %in% names(res_boot)) {
    res_boot <- add_survival_scale(res_boot, alpha)
  } else if (!"fragility_matched" %in% names(res_boot)) {
    res_boot <- add_survival_scale(res_boot, alpha)
  }
  
  q1 <- function(x) stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)
  q3 <- function(x) stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  
  out <- res_boot |>
    dplyr::filter(!is.na(delta_star_matched)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_iter = dplyr::n(),
      
      ## delta* matched (scala log-hazard)
      delta_star_matched_med = stats::median(delta_star_matched, na.rm = TRUE),
      delta_star_matched_q1  = q1(delta_star_matched),
      delta_star_matched_q3  = q3(delta_star_matched),
      abs_delta_star_matched_med = stats::median(abs(delta_star_matched), na.rm = TRUE),
      
      ## delta* registry (naive)
      delta_star_registry_med = stats::median(delta_star_registry, na.rm = TRUE),
      delta_star_registry_q1  = q1(delta_star_registry),
      delta_star_registry_q3  = q3(delta_star_registry),
      abs_delta_star_registry_med = stats::median(abs(delta_star_registry), na.rm = TRUE),
      
      ## tolleranza sulla scala sopravvivenza (S_tip - S_bench)
      d_star_matched_med = stats::median(d_star_matched, na.rm = TRUE),
      S_tip_matched_med  = stats::median(S_tip_matched,  na.rm = TRUE),
      
      ## rapporto di fragilita' (diagnostica di simulazione)
      fragility_matched_med  = stats::median(fragility_matched, na.rm = TRUE),
      frac_fragile_matched   = mean(fragility_matched > 1, na.rm = TRUE),
      frac_fragile_registry  = mean(fragility_registry > 1, na.rm = TRUE),
      
      ## contesto decisionale: frazione go implicata da delta*
      frac_go_matched  = mean(delta_star_matched  < 0, na.rm = TRUE),
      frac_go_registry = mean(delta_star_registry < 0, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  if (return_flextable) flextable::flextable(out) else out
}
