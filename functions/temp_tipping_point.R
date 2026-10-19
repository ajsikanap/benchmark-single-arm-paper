# tipping_point_tables.R
# ---------------------------------------------------------------------------
# Funzioni per le due tabelle del tipping point:
#   (1) tabella di DESIGN: theta*, delta*, S_tip(tau), m_tip su una griglia di
#       O plausibili, a E pianificato (nessun dato osservato richiesto);
#   (2) tabella dei RISULTATI di simulazione: distribuzione di delta* per
#       scenario x benchmark x N x effect_setting, confrontata con lo shift
#       vero tra benchmark e controfattuale (disponibile solo in simulazione).
#
# Notazione del paper (sez. 2.4-2.5):
#   L = (O - E) / sqrt(E), 
#   E = sum_i Lambda0(Y_i)
#   L(theta) = (O - theta*E) / sqrt(theta*E), 
#   L(theta*) = -z_alpha
#   theta* = (z_alpha + sqrt(z_alpha^2 + 4*O))^2 / (4*E) 
#   delta* = log(theta*)
#   GO  <=> theta* < 1 (delta* < 0)
#   S_tip(tau) = S0(tau)^theta*,  
#   m_tip = S0^{-1}(0.5^(1/theta*))
#   con S0^{-1}(p) = inf{t : S0(t) <= p}
#
# Non sostituisce tipping_point.R: le funzioni hanno nomi diversi.
# ---------------------------------------------------------------------------


# --- 0. Blocchi di base ------------------------------------------------------

tp_theta_star <- function(O, E, alpha = 0.10) {
  z <- stats::qnorm(1 - alpha)
  (z + sqrt(z^2 + 4 * O))^2 / (4 * E)
}

# Inversa generalizzata della survival del benchmark.
#   S0 = funzione continua t -> S0(t)          (es. esponenziale del protocollo)
#   S0 = data.frame(time, surv) a gradini      (es. Kaplan-Meier)

tp_inverse_surv <- function(p, S0, upper = 500) {
  if (is.function(S0)) {
    vapply(p, function(pp) {
      if (S0(upper) > pp) return(NA_real_)          # curva non scende a pp
      stats::uniroot(function(t) S0(t) - pp, c(0, upper))$root
    }, numeric(1))
  } else {
    vapply(p, function(pp) {
      idx <- which(S0$surv <= pp)
      if (length(idx) == 0) NA_real_ else min(S0$time[idx])
    }, numeric(1))
  }
}


tp_median_tip <- function(theta, S0) tp_inverse_surv(0.5^(1 / theta), S0)


# --- 1. Tabella di DESIGN ----------------------------------------------------

# E pianificato sotto il benchmark esponenziale del protocollo.
# Sotto H0 E[sum_i Lambda0(Y_i)] = E[numero di eventi] (proprieta' di
# martingala), quindi E_plan = n * P(T <= C), con T ~ Exp(log2/median0) e
# censura C ~ Uniform(C_low, C_upp) come nella simulazione.

tp_planned_E <- function(n, median0, C_low, C_upp) {
  lambda <- log(2) / median0
  p_event <- if (C_upp == C_low) {
    1 - exp(-lambda * C_upp)
  } else {
    1 - (exp(-lambda * C_low) - exp(-lambda * C_upp)) /
      (lambda * (C_upp - C_low))
  }
  n * p_event
}

tp_design_table <- function(E_plan,
                            O_grid  = NULL,
                            alpha   = 0.10,
                            tau     = 12,
                            median0 = 7.5,
                            S0      = NULL) {   # default: esponenziale del protocollo
  if (is.null(S0)) S0 <- function(t) exp(-log(2) / median0 * t)
  S0_tau <- if (is.function(S0)) S0(tau) else
    S0$surv[max(which(S0$time <= tau))]
  if (is.null(O_grid)) O_grid <- 0:ceiling(E_plan + 2 * sqrt(E_plan))
  
  z     <- stats::qnorm(1 - alpha)
  theta <- tp_theta_star(O_grid, E_plan, alpha)
  L     <- (O_grid - E_plan) / sqrt(E_plan)
  
  out <- data.frame(
    O           = O_grid,
    L           = L,
    decision    = ifelse(L < -z, "Go", "No-go"),
    theta_star  = theta,
    delta_star  = log(theta),
    S_tip_tau   = S0_tau^theta,
    median_tip  = tp_median_tip(theta, S0)
  )
  attr(out, "E_plan")  <- E_plan
  attr(out, "O_crit")  <- max(O_grid[L < -z])   # O massimo che porta a GO
  attr(out, "S0_tau")  <- S0_tau
  out
}

