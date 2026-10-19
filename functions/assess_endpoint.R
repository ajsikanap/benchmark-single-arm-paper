# function that 

assess_endpoint <- function(trial, matched, registry, tau = 12) {
  
  
  fit_trial_z1 <- survfit(Surv(time, event) ~ 1, data = trial)
  fit_trial_z0 <- survfit(Surv(time_z0, event_z0) ~ 1, data = trial)
  fit_matched  <- survfit(Surv(time, event) ~ 1, data = subset(matched, source == 0))
  fit_registry <- survfit(Surv(time, event) ~ 1, data = registry)
  
  log_rank_t<-survdiff(Surv(time,event)~source, data=matched)
  log_rank_pvalue=log_rank_t$pvalue
  
  
  data.frame(
    trial_z1 = get_surv_12(fit_trial_z1, time = tau),
    trial_true_z1 = mean(trial$T_true_z1 > tau, na.rm = TRUE),
    
    trial_z0 = get_surv_12(fit_trial_z0, time = tau),
    trial_true_z0 = mean(trial$T_true_z0 > tau, na.rm = TRUE),
    
    matched_z0  = get_surv_12(fit_matched,  time = tau),
    matched_true_z0 = mean(subset(matched, source == 0)$T_true_z0>tau, na.rm=TRUE),
    
    registry = get_surv_12(fit_registry, time = tau),
    registry_true = mean(registry$T_true_z0>tau, na.rm=TRUE),
    
    rmst_trial_z1 = get_rmst_from_km(fit_trial_z1, tau = tau),
    rmst_trial_true_z1=mean(pmin(trial$T_true_z1, tau), na.rm = TRUE),
    
    rmst_trial_z0 = get_rmst_from_km(fit_trial_z0, tau = tau),
    rmst_trial_true_z0=mean(pmin(trial$T_true_z0, tau), na.rm = TRUE),
    
    rmst_matched_z0  = get_rmst_from_km(fit_matched,  tau = tau),
    rmst_matched_true_z0 = mean(pmin(subset(matched, source == 0)$T_true_z0,tau), na.rm=TRUE),
    
    rmst_registry = get_rmst_from_km(fit_registry, tau = tau),
    rmst_registry_true=mean(pmin(registry$T_true_z0, tau), na.rm = TRUE),
    
    lr_two_sample_pvalue=log_rank_pvalue
  )
}


