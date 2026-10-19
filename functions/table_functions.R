# =============================================================================
# table_functions.R
# Functions to build Tables 3-6b for the paper.
# Source this file in 06_paper_outputs.Rmd after loading results_raw.RData
# and sourcing load_functions.R (che a sua volta sorge summarise_main_results.R).
#
# AGGIORNAMENTO: summarise_main_results() e' stata estesa per produrre in un
# unico data frame bias/RMSE, decomposizione del bias, livelli assoluti
# (latenti e stimati) e operating characteristics. Le funzioni sotto sono
# state riallineate ai nomi di colonna che quella funzione produce ora.
# Le vecchie versioni di build_table3/4/5 (pensate per balance_all/all_results/
# bias_decomposition_all separati) sono state rimosse perche' non corrispondono
# piu' alla pipeline attuale; il codice e' recuperabile dalla history di git
# se serve confrontarlo.
# =============================================================================

library(dplyr)
library(tidyr)
library(knitr)
library(flextable)
library(ggplot2)

# ── Shared helpers ───────────────────────────────────────────────────────────

# Scenario label ordering for display
SCENARIO_ORDER <- c("S0", "S1_all", "S2", "S3", "S5", "S4")

SCENARIO_LABELS <- c(
  S0     = "S0 (oracle)",
  S1_all = "S1",
  S2     = "S2",
  S3     = "S3A",
  S5     = "S3B",
  S4     = "S4"
)

METHOD_LABELS <- c(prop = "Propensity", prog = "Prognostic")

fmt <- function(x, d = 3) round(x, d)

# Format mean (SD) as a single string
fmt_mean_sd <- function(mean_col, sd_col, d = 3) {
  paste0(fmt(mean_col, d), " (", fmt(sd_col, d), ")")
}

# Formatta "stima (mcse)" in una cella unica.
# d_est / d_mcse: decimali per stima e per MCSE (il MCSE spesso vuole 1 cifra in piu').
fmt_est_mcse <- function(est, mcse, d_est = 3, d_mcse = 4) {
  ifelse(
    is.na(est), "\u2014",
    paste0(formatC(est,  format = "f", digits = d_est),
           " (", formatC(mcse, format = "f", digits = d_mcse), ")")
  )
}

# Rename scenario_block using SCENARIO_LABELS
relabel_scenario <- function(df) {
  df |>
    dplyr::mutate(
      scenario_block = factor(
        dplyr::recode(scenario_block, !!!SCENARIO_LABELS),
        levels = unname(SCENARIO_LABELS)
      ),
      registry_n = factor(registry_n, levels = c(75,100,125,175, 250, 500, 1000)),
      method     = dplyr::recode(method, !!!METHOD_LABELS)
    )
}

# Standard flextable formatting
paper_flextable <- function(df, caption_text) {
  df |>
    flextable::flextable() |>
    flextable::set_caption(caption_text) |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::bg(bg = "#D0E4F0", part = "header") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1, align = "left", part = "all") |>
    flextable::border_outer(border = officer::fp_border(color = "black", width = 1)) |>
    flextable::border_inner_h(border = officer::fp_border(color = "#BBBBBB", width = 0.5)) |>
    flextable::autofit()
}

# =============================================================================
# TABLE 3 - Bias e RMSE dell'effetto di trattamento (formato "stima (MCSE)")
# Estimand primario: DeltaS(tau).  Co-primario: DeltaRMST(tau).
# Source: summarise_main_results(res_boot_alt, group_vars = c("scenario","analysis_label","registry_n","method"))
#         con "scenario" rinominato in "scenario_block" (vedi to_block() nel .Rmd).
#
# NB sui nomi: la stima puntuale del bias si chiama total_bias_deltaS_matched /
# total_bias_deltaS_naive (non bias_deltaS_matched / bias_deltaS_registry), mentre
# il relativo MCSE ha mantenuto il vecchio suffisso (bias_deltaS_matched_mcse /
# bias_deltaS_registry_mcse). E' un'asimmetria ereditata da summarise_main_results();
# se la cambi la' di nuovo, aggiorna anche qui.
# =============================================================================

build_table3 <- function(main_summary, scale = c("deltaS", "deltaRMST")) {

  scale <- match.arg(scale)

  if (scale == "deltaS") {
    cap <- paste0(
      "Table 3. Bias and RMSE of the estimated treatment effect on 12-month ",
      "survival difference (\u0394S), by scenario, registry size, and method. ",
      "Matched: score-matched benchmark; Naive: full-registry benchmark. ",
      "Monte Carlo standard errors in parentheses. n = valid replicates ",
      "(matching successful)."
    )
    tab <- main_summary |>
      dplyr::transmute(
        scenario_block, registry_n, method, n_valid,
        `Bias matched`  = fmt_est_mcse(total_bias_deltaS_matched, bias_deltaS_matched_mcse),
        `RMSE matched`  = fmt_est_mcse(rmse_deltaS_matched,       rmse_deltaS_matched_mcse),
        `Bias naive`    = fmt_est_mcse(total_bias_deltaS_naive,   bias_deltaS_registry_mcse),
        `RMSE naive`    = fmt_est_mcse(rmse_deltaS_registry,      rmse_deltaS_registry_mcse)
      )
  } else {
    cap <- paste0(
      "Table 3 (co-primary). Bias and RMSE of the estimated treatment effect on ",
      "restricted mean survival time (\u0394RMST, \u03c4 = 12), by scenario, registry ",
      "size, and method. Matched vs Naive benchmark. Monte Carlo standard errors ",
      "in parentheses."
    )
    tab <- main_summary |>
      dplyr::transmute(
        scenario_block, registry_n, method, n_valid,
        `Bias matched`  = fmt_est_mcse(total_bias_deltaRMST_matched, bias_deltaRMST_matched_mcse,  d_est = 3, d_mcse = 3),
        `RMSE matched`  = fmt_est_mcse(rmse_deltaRMST_matched,       rmse_deltaRMST_matched_mcse,  d_est = 3, d_mcse = 3),
        `Bias naive`    = fmt_est_mcse(total_bias_deltaRMST_naive,   bias_deltaRMST_registry_mcse, d_est = 3, d_mcse = 3),
        `RMSE naive`    = fmt_est_mcse(rmse_deltaRMST_registry,      rmse_deltaRMST_registry_mcse, d_est = 3, d_mcse = 3)
      )
  }

  tab |>
    relabel_scenario() |>
    dplyr::arrange(scenario_block, registry_n, method) |>
    dplyr::select(
      Scenario     = scenario_block,
      `Registry N` = registry_n,
      Method       = method,
      n            = n_valid,
      `Bias matched`, `RMSE matched`, `Bias naive`, `RMSE naive`
    ) |>
    paper_flextable(cap)
}

# =============================================================================
# build_table_main() - Tabella estesa: livelli veri (latenti) + differenza
# stimata + bias/RMSE, naive e matched fianco a fianco. Pensata come tabella
# "quadro completo" da supplementary, complementare a build_table3 (solo
# bias/RMSE) e build_table6a (solo livelli/medie assolute, senza incertezza).
#
# NB di disegno: Trial Z=1 / Trial Z=0 / Naive / Matched sono i livelli VERI
# (latenti, dalla superpopulation), non le stime KM osservate - servono da
# contesto descrittivo. "Estimated Delta" e "Total Bias"/"RMSE" sono invece
# sulla scala osservata/stimata. Se preferisci i livelli osservati (KM) al
# posto di quelli latenti, sostituisci mean_registry_true/mean_matched_true_z0
# con mean_registry_est/mean_matched_est (e gli equivalenti RMST).
# =============================================================================

build_table_main <- function(main_summary, scale = c("deltaS", "deltaRMST"), d = 3) {

  scale <- match.arg(scale)
  fmt_d <- function(x) formatC(round(x, d), format = "f", digits = d)

  if (scale == "deltaS") {
    cap <- ""
    tab <- main_summary |>
      dplyr::transmute(
        scenario_block, registry_n, method, n_valid,
        trial_z1          = mean_trial_true_z1,
        trial_z0_true      = mean_trial_true_z0,
        naive_level         = mean_registry_true,
        matched_level       = mean_matched_true_z0,
        true_delta          = true_deltaS,
        est_delta_naive     = mean_deltaS_naive_obs,
        bias_naive          = fmt_est_mcse(total_bias_deltaS_naive,  bias_deltaS_registry_mcse),
        rmse_naive          = fmt_est_mcse(rmse_deltaS_registry,     rmse_deltaS_registry_mcse),
        est_delta_matched   = mean_deltaS_matched,
        bias_matched        = fmt_est_mcse(total_bias_deltaS_matched, bias_deltaS_matched_mcse),
        rmse_matched        = fmt_est_mcse(rmse_deltaS_matched,       rmse_deltaS_matched_mcse)
      )
  } else {
    cap <- ""
    tab <- main_summary |>
      dplyr::transmute(
        scenario_block, registry_n, method, n_valid,
        trial_z1          = mean_rmst_trial_true_z1,
        trial_z0_true      = mean_rmst_trial_true_z0,
        naive_level         = mean_rmst_registry_true,
        matched_level       = mean_rmst_matched_true_z0,
        true_delta          = true_deltaRMST,
        est_delta_naive     = mean_deltaRMST_naive_obs,
        bias_naive          = fmt_est_mcse(total_bias_deltaRMST_naive,  bias_deltaRMST_registry_mcse, d_est = 3, d_mcse = 3),
        rmse_naive          = fmt_est_mcse(rmse_deltaRMST_registry,     rmse_deltaRMST_registry_mcse, d_est = 3, d_mcse = 3),
        est_delta_matched   = mean_deltaRMST_matched,
        bias_matched        = fmt_est_mcse(total_bias_deltaRMST_matched, bias_deltaRMST_matched_mcse, d_est = 3, d_mcse = 3),
        rmse_matched        = fmt_est_mcse(rmse_deltaRMST_matched,       rmse_deltaRMST_matched_mcse, d_est = 3, d_mcse = 3)
      )
  }

  # Solo le colonne realmente numeriche vengono arrotondate qui; le colonne
  # bias_*/rmse_* sono GIA' stringhe formattate da fmt_est_mcse() e non vanno
  # passate a fmt_d() (altrimenti round() su un vettore character va in errore).
  num_cols <- c("trial_z1", "trial_z0_true", "naive_level", "matched_level",
                "true_delta", "est_delta_naive", "est_delta_matched")

  tab |>
    relabel_scenario() |>
    dplyr::arrange(scenario_block, registry_n, method) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(num_cols), fmt_d)) |>
    dplyr::select(
      Scenario             = scenario_block,
      `Registry N`         = registry_n,
      Method               = method,
      n                    = n_valid,
      `Trial Z=1`          = trial_z1,
      `Trial Z=0 (true)`   = trial_z0_true,
      `Naive`              = naive_level,
      `Matched`            = matched_level,
      `True Delta`             = true_delta,
      `Estimated  Delta Naive`  = est_delta_naive,
      `Total Bias Naive`   = bias_naive,
      `RMSE Naive`         = rmse_naive,
      `Estimate Delta Matched` = est_delta_matched,
      `Total Bias Matched` = bias_matched,
      `RMSE Matched`       = rmse_matched
    ) |>
    paper_flextable(cap)
}

# =============================================================================
# build_table4() - Decomposizione del bias dell'effetto di trattamento.
#
# Estimand:  scale = "deltaS"  (primario) | "deltaRMST" (co-primario)
# Benchmark: bench = "matched" (il metodo)| "naive"     (registry, confronto)
#
# Lettura per riga (matched):
#   Total bias = Trial est.err - Matched est.err - Matching err + Sampling err
# (naive):
#   Total bias = Trial est.err - Registry est.err - Selection err + Sampling err
#
# "Check" e' il residuo algebrico (Total bias - ricostruito dalla somma dei
# componenti) e va a ~0 quando la decomposizione e' esatta. Va calcolato qui
# da total_bias_* e reconstructed_bias_* (summarise_main_results() non
# espone piu' una colonna check_* pre-calcolata).
# =============================================================================

build_table4 <- function(decomp_summary,
                          scale = c("deltaS", "deltaRMST"),
                          bench = c("matched", "naive"),
                          d = 3) {

  scale <- match.arg(scale)
  bench <- match.arg(bench)

  fmt_d <- function(x) formatC(round(x, d), format = "f", digits = d)

  # ── Selezione colonne in base a scale x bench ──────────────────────────────
  if (scale == "deltaS" && bench == "matched") {
    tab <- decomp_summary |>
      dplyr::transmute(
        scenario, analysis_label, registry_n, method,
        `Mean treatment difference` = mean_deltaS_matched_obs,
        `True treatment difference` = true_deltaS,
        `Total bias`          = total_bias_deltaS_matched,
        `Trial est. error`    = trial_est_error_S,
        `-Control est. error` = -matched_control_est_error_S,
        `-Matching error`     = -matching_error_S,
        `Sampling error`      = sampling_error_deltaS,
        `Check`               = total_bias_deltaS_matched - reconstructed_bias_deltaS_matched
      )
    cap_scale <- "12-month survival difference (\u0394S)"; cap_bench <- "matched-control"

  } else if (scale == "deltaS" && bench == "naive") {
    tab <- decomp_summary |>
      dplyr::transmute(
        scenario, analysis_label, registry_n, method,
        `Mean treatment difference` = mean_deltaS_naive_obs,
        `True treatment difference` = true_deltaS,
        `Total bias`           = total_bias_deltaS_naive,
        `Trial est. error`     = trial_est_error_S,
        `-Registry est. error` = -registry_est_error_S,
        `-Selection error`     = -registry_selection_error_S,
        `Sampling error`       = sampling_error_deltaS,
        `Check`                = total_bias_deltaS_naive - reconstructed_bias_deltaS_naive
      )
    cap_scale <- "12-month survival difference (\u0394S)"; cap_bench <- "full-registry (naive)"

  } else if (scale == "deltaRMST" && bench == "matched") {
    tab <- decomp_summary |>
      dplyr::transmute(
        scenario, analysis_label, registry_n, method,
        `Mean treatment difference` = mean_deltaRMST_matched_obs,
        `True treatment difference` = true_deltaRMST,
        `Total bias`          = total_bias_deltaRMST_matched,
        `Trial est. error`    = trial_est_error_RMST,
        `-Control est. error` = -matched_control_est_error_RMST,
        `-Matching error`     = -matching_error_RMST,
        `Sampling error`      = sampling_error_deltaRMST,
        `Check`               = total_bias_deltaRMST_matched - reconstructed_bias_deltaRMST_matched
      )
    cap_scale <- "restricted mean survival time (\u0394RMST)"; cap_bench <- "matched-control"

  } else {  # deltaRMST naive
    tab <- decomp_summary |>
      dplyr::transmute(
        scenario, analysis_label, registry_n, method,
        `Mean treatment difference` = mean_deltaRMST_naive_obs,
        `True treatment difference` = true_deltaRMST,
        `Total bias`           = total_bias_deltaRMST_naive,
        `Trial est. error`     = trial_est_error_RMST,
        `-Registry est. error` = -registry_est_error_RMST,
        `-Selection error`     = -registry_selection_error_RMST,
        `Sampling error`       = sampling_error_deltaRMST,
        `Check`                = total_bias_deltaRMST_naive - reconstructed_bias_deltaRMST_naive
      )
    cap_scale <- "restricted mean survival time (\u0394RMST)"; cap_bench <- "full-registry (naive)"
  }

  # ── Ordina/rietichetta, formatta i numeri, rinomina colonne chiave ─────────
  num_cols <- setdiff(names(tab), c("scenario", "analysis_label", "registry_n", "method"))

  tab |>
    dplyr::rename(scenario_block = scenario) |>
    relabel_scenario() |>
    dplyr::arrange(scenario_block, registry_n, method) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(num_cols), fmt_d)) |>
    dplyr::rename(Scenario = dplyr::any_of(c("scenario_block", "scenario"))) |>
    dplyr::select(
      Scenario,
      `Registry N` = registry_n,
      Method       = method,
      dplyr::all_of(num_cols)
    ) |>
    paper_flextable(
      paste0(
        "Table 4. Decomposition of the bias in the estimated treatment effect (with sign integrated in the value) on ",
        cap_scale, ", using the ", cap_bench, " benchmark, by scenario and registry ",
        "size. The total bias decomposes additively into estimation error in the ",
        "trial arm, estimation error in the control benchmark, the ",
        if (bench == "matched") "matching error (residual confounding after matching), "
        else "registry selection error (non-transportability of the historical control), ",
        "and the Monte Carlo sampling error. 'Check' is the algebraic residual ",
        "(total minus reconstructed) and is \u2248 0 when the decomposition is exact. ",
        "Values are means across replicates with successful matching."
      )
    )
}

# =============================================================================
# TABLE 5 - Operating characteristics
# Tasso di rigetto (one-sided, alpha = 0.10) per scenario x registry x metodo.
# Da res_boot_NULL -> Type I error ; da res_boot_ALT -> Power.
# Colonne: one-sample matched, one-sample naive, two-sample matched, two-sample naive.
# Source: summarise_main_results() su null e su alt, poi passati qui insieme.
# Nessuna modifica necessaria qui: i nomi reject_* corrispondono gia' a quelli
# prodotti da summarise_main_results().
# =============================================================================

build_table5 <- function(null_summary, alt_summary) {

  pick <- function(df, label) {
    df |>
      dplyr::transmute(
        scenario_block, registry_n, method,
        !!paste(label, "One-sample (matched)") := fmt_est_mcse(reject_one_sample_matched, reject_one_sample_matched_mcse),
        !!paste(label, "One-sample (naive)") := fmt_est_mcse(reject_one_sample_registry, reject_one_sample_registry_mcse),
        !!paste(label, "Two-sample (matched)") := fmt_est_mcse(reject_two_sample_matched, reject_two_sample_matched_mcse),
        !!paste(label, "Two-sample (naive)")   := fmt_est_mcse(reject_two_sample_naive,   reject_two_sample_naive_mcse)
      )
  }

  dplyr::left_join(
    pick(null_summary, "Type I error"),
    pick(alt_summary,  "Power"),
    by = c("scenario_block", "registry_n", "method")
  ) |>
    relabel_scenario() |>
    dplyr::arrange(scenario_block, registry_n, method) |>
    dplyr::select(
      Scenario = scenario_block,
      `Registry N` = registry_n,
      Method       = method,
      `Type I error One-sample (matched)`, `Type I error One-sample (naive)`,
      `Type I error Two-sample (matched)`, `Type I error Two-sample (naive)`,
      `Power One-sample (matched)`, `Power One-sample (naive)`,
      `Power Two-sample (matched)`, `Power Two-sample (naive)`
    ) |>
    paper_flextable(
      paste0(
        "Table 5. Operating characteristics: one-sided rejection rate at ",
        "\u03b1 = 0.10. Type I error from the null setting (no treatment effect), ",
        "power from the alternative. One-sample log-rank treats the benchmark as ",
        "fixed; two-sample compares against matched controls (matched) or the full ",
        "registry (naive). Monte Carlo standard errors in parentheses. Type I error ",
        "inflation above 0.10 signals a non-transportable benchmark."
      )
    )
}

# =============================================================================
# TABLE 6a - Mean benchmark estimates (Delta S and Delta RMST)
# Source: decision_summary_all
#
# NB: tutte le colonne richieste qui (mean_deltaS_oracle/registry/matched,
# mean_trial_z1_est/mean_trial_z0_est/mean_matched_est/mean_registry_est,
# mean_deltaRMST_*, mean_rmst_*) sono prodotte direttamente da
# summarise_main_results(). Se "decision_summary_all" non e' ancora costruito
# nel .Rmd, in pratica puoi chiamare build_table6a(table_main_alt) cosi' com'e'
# (non serve relabel_scenario/to_block: la funzione rietichetta scenario al
# suo interno con la sua propria mutate()).
# =============================================================================

build_table6a <- function(decision_summary_all) {

  decision_summary_all |>
    dplyr::mutate(
      scenario = ifelse(scenario == "S0", "S0 (oracle)", scenario),
      scenario = ifelse(scenario == "S3", "S3A", scenario),
      scenario = ifelse(scenario == "S5", "S3B", scenario),
      scenario   = factor(scenario, levels = unname(SCENARIO_LABELS)),
      registry_n = factor(registry_n, levels = c(125, 250, 500, 1000)),
      method     = dplyr::recode(method, !!!METHOD_LABELS)
    ) |>
    dplyr::arrange(scenario, registry_n, method) |>
    dplyr::select(
      Scenario              = scenario,
      `Registry N`          = registry_n,
      Method                = method,
      # Survival difference
      ` Delta S oracle`       = mean_deltaS_oracle,
      ` Delta S registry`     = mean_deltaS_registry,
      ` Delta S matched`      = mean_deltaS_matched,

      # Survival levels
      `S(t) trial Z=1`     = mean_trial_z1_est,
      `S(t) trial Z=0`     = mean_trial_z0_est,
      `S(t) matched`   = mean_matched_est,
      `S(t) registry`  = mean_registry_est,

      # rmst difference
      ` Delta RMST oracle`       = mean_deltaRMST_oracle,
      ` Delta RMST registry`     = mean_deltaRMST_registry,
      ` Delta RMST matched`      = mean_deltaRMST_matched,

      # RMST
      `RMST trial Z=1`     = mean_rmst_trial_z1,
      `RMST trial Z=0`     = mean_rmst_trial_z0,
      `RMST matched`   = mean_rmst_matched,
      `RMST registry`  = mean_rmst_registry
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) fmt(x, 3))) |>
    paper_flextable(
      "Table 6a. Mean estimated treatment effects under each benchmark definition,
       by scenario and registry size. Delta S: survival difference at \u03c4 = 12 months;
       RMST: restricted mean survival time up to \u03c4 = 12 months.
       Oracle benchmark uses the true sampled trial counterfactual (not available in practice).
       Values are means across simulation replicates."
    )
}

# =============================================================================
# TABLE 6b - Operating characteristics (Type I error and Power)
# Source: decision_summary_all
#
# NB: a differenza di 6a, QUESTA tabella si aspetta colonne con la
# convenzione di nomi "typeI_one_sample_lr_registry_1s" /
# "power_two_sample_lr_matched_1s" ecc., che NON corrispondono ai nomi
# reject_one_sample_matched / reject_two_sample_naive prodotti oggi da
# summarise_main_results() (vedi build_table5 e table_oc nel .Rmd, che usano
# la convenzione attuale). E' probabile un residuo di una versione precedente
# della pipeline. Prima di usarla, va riscritta sul modello di build_table5/
# table_oc, oppure va costruito un "decision_summary_all" con questi nomi
# specifici - non l'ho rifatta perche' senza vedere quale dei due percorsi
# preferisci rischierei di introdurre una terza convenzione di naming.
# =============================================================================

build_table6b <- function(decision_summary_all) {

  decision_summary_all |>
    dplyr::mutate(
      scenario = ifelse(scenario == "S0", "S0 (oracle)", scenario),
      scenario = ifelse(scenario == "S3", "S3A", scenario),
      scenario = ifelse(scenario == "S5", "S3B", scenario),
      scenario   = factor(scenario, levels = unname(SCENARIO_LABELS)),
      registry_n = factor(registry_n, levels = c(125, 250, 500, 1000)),
      method     = dplyr::recode(method, !!!METHOD_LABELS)
    ) |>
    dplyr::arrange(scenario, registry_n, method) |>
    dplyr::select(
      Scenario       = scenario,
      `Registry N`   = registry_n,
      Method         = method,
      # One-sample log-rank
      `Type I\n1s-LR registry`  = typeI_one_sample_lr_registry_1s,
      `Type I\n1s-LR matched`   = typeI_one_sample_lr_matched_1s,
      `Power\n1s-LR registry`   = power_one_sample_lr_registry_1s,
      `Power\n1s-LR matched`    = power_one_sample_lr_matched_1s,
      # Two-sample log-rank
      `Type I\n2s-LR registry`  = typeI_two_sample_lr_registry_1s,
      `Type I\n2s-LR matched`   = typeI_two_sample_lr_matched_1s,
      `Power\n2s-LR registry`   = power_two_sample_lr_registry_1s,
      `Power\n2s-LR matched`    = power_two_sample_lr_matched_1s
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) fmt(x, 3))) |>
    paper_flextable(
      "Table 6b. Operating characteristics of one-sample (1s-LR) and two-sample (2s-LR)
       log-rank tests under each benchmark definition, by scenario and registry size.
       Type I error: probability of rejecting H under the null scenario.
       Power: probability of rejecting H under the alternative scenario.
       Nominal \u03b1 = 0.10 (one-sided); target power = 0.80.
       Registry: na\u00efve full registry benchmark; matched: propensity or prognostic matched."
    )
}

# =============================================================================
# USAGE - vedi 06_paper_outputs.Rmd per l'ordine effettivo delle chiamate.
#
# to_block <- function(df) dplyr::rename(df, scenario_block = scenario)
#
# build_table3(to_block(table_main_alt), scale = "deltaS")
# build_table3(to_block(table_main_alt), scale = "deltaRMST")
# build_table_main(to_block(table_main_alt), scale = "deltaS")
# build_table_main(to_block(table_main_alt), scale = "deltaRMST")
#
# build_table4(table_main_alt, scale = "deltaS",    bench = "matched")
# build_table4(table_main_alt, scale = "deltaS",    bench = "naive")
# build_table4(table_main_alt, scale = "deltaRMST", bench = "matched")
# build_table4(table_main_alt, scale = "deltaRMST", bench = "naive")
#
# build_table5(to_block(table_main_null), to_block(table_main_alt))
#
# build_table6a(table_main_alt)   # gia' compatibile, vedi nota sopra
# build_table6b(decision_summary_all)  # da rivedere, vedi nota sopra
# =============================================================================
