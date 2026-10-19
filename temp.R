# =============================================================================
# check_refactor.R  —  temporary check, NOT part of the repository
#
# Level 1: old 04 script vs cleaned 04 script, same N_BOOT, same environment.
#          Expected: identical results (isolates the effect of the cleanup).
# Level 2: cleaned 04 (first N_BOOT iterations) vs the paper run
#          results_raw_10000.RData restricted to the same iterations.
#          Iteration i depends only on seeds[i], so it must not depend on
#          N_BOOT. Differences here come from functions/packages that changed
#          since the paper run, not from the cleanup.
# =============================================================================

key_cols <- c("scenario", "registry_n", "method", "effect_setting", "iter")

load_env <- function(file) {
  e <- new.env()
  load(file, envir = e)
  e
}

prep <- function(df) {
  df <- as.data.frame(df)
  if (anyDuplicated(df[key_cols])) stop("Duplicated keys: key_cols do not identify rows.")
  df <- df[do.call(order, unname(as.list(df[key_cols]))), , drop = FALSE]
  rownames(df) <- NULL
  df
}

# Number of positions where two vectors differ (NA == NA counts as equal)
n_diff <- function(x, y) {
  both_na <- is.na(x) & is.na(y)
  same    <- both_na | (!is.na(x) & !is.na(y) & x == y)
  sum(!same)
}

compare_df <- function(a, b, label) {
  a <- prep(a); b <- prep(b)
  cat("\n──", label, "──\n")
  
  only_a <- setdiff(names(a), names(b))
  only_b <- setdiff(names(b), names(a))
  if (length(only_a)) cat("Columns only in first :", toString(only_a), "\n")
  if (length(only_b)) cat("Columns only in second:", toString(only_b), "\n")
  
  if (nrow(a) != nrow(b) || !identical(a[key_cols], b[key_cols])) {
    cat("Row sets differ:", nrow(a), "vs", nrow(b), "rows\n")
    ka <- do.call(paste, a[key_cols]); kb <- do.call(paste, b[key_cols])
    cat("  keys only in first :", sum(!ka %in% kb), "\n")
    cat("  keys only in second:", sum(!kb %in% ka), "\n")
    return(invisible(FALSE))
  }
  
  common   <- intersect(names(a), names(b))
  not_iden <- common[!vapply(common, \(v) identical(a[[v]], b[[v]]), logical(1))]
  
  if (!length(not_iden)) {
    cat("IDENTICAL:", nrow(a), "rows,", length(common), "common columns\n")
    return(invisible(TRUE))
  }
  
  out <- data.frame(
    column       = not_iden,
    class_first  = vapply(not_iden, \(v) class(a[[v]])[1], ""),
    class_second = vapply(not_iden, \(v) class(b[[v]])[1], ""),
    n_diff       = vapply(not_iden, \(v) n_diff(a[[v]], b[[v]]), 0L),
    max_abs_diff = vapply(not_iden, \(v) {
      if (is.numeric(a[[v]]) && is.numeric(b[[v]]))
        suppressWarnings(max(abs(a[[v]] - b[[v]]), na.rm = TRUE)) else NA_real_
    }, 0),
    row.names = NULL
  )
  # n_diff == 0 with identical() FALSE means only type/attributes differ
  print(out)
  invisible(out)
}

# ── Level 1: old vs cleaned script ────────────────────────────────────────────
new <- load_env(here::here("results/results_raw_10.RData"))
old <- load_env(("C:/Users/ajsi.kanapari/The Institute of Cancer Research/ICR-CTSU Early Phase and Adaptive Trials - DCS - Ajsi_external_control/historicalcontrols/results/results_raw_10.RData"))
res_boot_null[1:10,]

compare_df(old$res_boot_alt,  new$res_boot_alt,  "L1 alternative: old vs new")
compare_df(old$res_boot_null, new$res_boot_null, "L1 null: old vs new")
cat("\ntruth_scenarios identical :", identical(old$truth_scenarios,  new$truth_scenarios), "\n")
cat("scenario_configs identical:", identical(old$scenario_configs, new$scenario_configs), "\n")
cat("crashed iterations (new)  :", nrow(new$res_boot_errors), "\n")



# ── Level 2: cleaned script vs paper run (same iterations) ────────────────────
ref <- load_env(("C:/Users/ajsi.kanapari/The Institute of Cancer Research/ICR-CTSU Early Phase and Adaptive Trials - DCS - Ajsi_external_control/historicalcontrols/results/results_raw_10000.RData"))
k   <- new$N_BOOT

compare_df(subset(ref$res_boot_alt,  iter <= k), new$res_boot_alt,  "L2 alternative: paper vs new")
compare_df(subset(ref$res_boot_null, iter <= k), new$res_boot_null, "L2 null: paper vs new")

subset(ref$res_boot_alt, iter<=k)[c(1:5),c(2,33:53)]
new$res_boot_alt[c(1:5),c(2,33:53)]

flags <- c("M","age","sex","resection","trt_type","WHOm2","WHO","multifocal",
           "pTERT","tumor_size","EGFR","g7l10","unneccessary")
for (s in c("S1", "S4")) {
  get_flags <- \(cfg) dplyr::filter(cfg$sc, parameters == "trial") |> dplyr::select(dplyr::any_of(flags))
  print(dplyr::bind_rows(paper = get_flags(ref$scenario_configs[[s]]),
                         new   = get_flags(new$scenario_configs[[s]]), .id = s))
}

mean(registry_generated$B)                        # prevalenza effettiva nei dati del paper
c(registry = max(registry_generated$time),        # censura amministrativa: deve coincidere
  trial    = max(trial_pool_scenarios$S1$time))   # se C_upp e C_upp_trial sono unificati

a <- prep(subset(ref$res_boot_alt, iter <= k)); b <- prep(new$res_boot_alt)
same_abs <- abs(abs(a$score_smd_before) - b$score_smd_before) < 1e-8
tapply(same_abs, a$method, mean)   # quota di righe uguali per metodo
packageVersion("MatchIt")          # da lanciare in entrambi i progetti
