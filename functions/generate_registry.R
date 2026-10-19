# use L only if requested
# L only affcts trial inclusion
# L is a latent selection variable
# 
generate_registry <- function(
    median_z0
    ,N
    ,lambda0=log(2)/median_z0
    ,C_low
    ,C_upp
    
    
    # prevalences of prognostic covariates
    ,p_M
    ,p_B 
    ,p_sexF
    ,p_resection_0_1_2_3
    ,p_trt_baseline_0_1_2
    ,p_WHOm2
    ,p_loc_superficial_0
    ,p_loc_superficial_1
    ,p_loc_superficial_2
    ,p_loc_superficial_3
    ,p_loc_superficial_4
    ,p_loc_superficial_5
    ,p_multifocal
    ,mean_tumor_size
    ,sd_tumor_size
    ,p_pTERT 
    ,p_EGFR
    ,p_g7l10
    ,mean_age
    ,sd_age
    
    
    
    # mean and sd of unmeasured and unneccesary
    ,Unmeas_mean 
    ,Unmeas_sd 
    ,Unnec_mean 
    ,Unnec_sd
    ,p_unnec
    ,p_unmeas
    
    # beta - prog effect on OS
    ,beta_unmeas_binary
    ,beta_treatment_OS 
    ,int_OS 
    ,beta_M_OS
    ,beta_B_OS 
    ,beta_age_OS 
    ,beta_WHOm2_OS
    ,beta_sex_OS
    ,beta_resection_OS
    ,beta_trt_baseline_OS
    ,beta_localization_OS
    ,beta_multifocal
    ,beta_pTERT_OS
    ,beta_tumor_size_OS
    ,beta_EGFR_OS
    ,beta_g7l10_OS
    ,beta_Unmeas_OS
    
    # correlation
    ,rho_B_TS  # desired correlation between M and B
    ,rho_U_TS # correlation between unobserved and tumor size
    ,delta_TS
    ,delta_EGFR 
    ,use_L = FALSE
    ,L_mean
    ,L_sd 
    ,alpha0 
    ,alpha_B 
    ,alpha_L 
    ,alpha_U
    ,rho_L_age
    
) {
  N <- as.integer(N)
  
  rho_B_TS <- as.numeric(rho_B_TS)[1]
  rho_U_TS <- as.numeric(rho_U_TS)[1]
  
  rho_B_TS <- max(min(rho_B_TS, 0.999), -0.999)
  rho_U_TS <- max(min(rho_U_TS, 0.999), -0.999)
  
  Sigma_B_TS <- matrix(
    c(1, rho_B_TS,
      rho_B_TS, 1),
    nrow = 2,
    byrow = TRUE
  )  
  
  
  latent <- MASS::mvrnorm(
    n = N,
    mu = c(0, 0),
    Sigma = Sigma_B_TS
  )
  
  TS_latent <- latent[, 1]
  B_latent  <- latent[, 2]
  
  tumor_size <- mean_tumor_size + sd_tumor_size * TS_latent #
  tumor_size <- pmax(tumor_size, 1)
  
  # for simplicity they move in the same direction
  B <- as.numeric(B_latent > qnorm(p_B)) 
  
  # Correlated latent variables for unmeasured binary component.
  Sigma_U_TS <- matrix(c(1, rho_U_TS, rho_U_TS, 1), nrow = 2)
  latent1 <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_U_TS)
  
  unmeas_binary <- as.numeric(latent1[, 2] < qnorm(p_unmeas))
  
  registry<-data.frame(
    id = 1:N,
    M = rbinom(N,1,p_M),
    B = B,
    age = rnorm(N, mean_age, sd_age),
    sex = rbinom(N, 1, p_sexF),
    WHOm2 = rbinom(N, 1, p_WHOm2),
    resection = sample(0:3, N, TRUE, p_resection_0_1_2_3),
    trt_type = sample(0:2, N, TRUE, p_trt_baseline_0_1_2),
    localization_0 = rbinom(N, 1, p_loc_superficial_0),
    localization_1 = rbinom(N, 1, p_loc_superficial_1),
    localization_2 = rbinom(N, 1, p_loc_superficial_2),
    localization_3 = rbinom(N, 1, p_loc_superficial_3),
    localization_4 = rbinom(N, 1, p_loc_superficial_4),
    localization_5 = rbinom(N, 1, p_loc_superficial_5),
    multifocal = rbinom(N, 1, p_multifocal),
    tumor_size = tumor_size,
    pTERT = rbinom(N, 1, p_pTERT),
    EGFR = rbinom(N, 1, p_EGFR),
    g7l10 = rbinom(N, 1, p_g7l10),
    unmeas_binary = unmeas_binary,
    unmeasured = rnorm(N, Unmeas_mean, Unmeas_sd),
    unneccessary = rbinom(N, 1, p_unnec)
  )
  
  
  
  ## OS
  C <- runif(N, C_low, C_upp)
  
  lp_OS_z0 <- with(registry,
                   int_OS +
                     beta_M_OS*M +
                     beta_B_OS*B +
                     beta_age_OS*(age - 60) +
                     beta_Unmeas_OS*unmeasured +
                     WHOm2*beta_WHOm2_OS +
                     sex*beta_sex_OS +
                     resection*beta_resection_OS +
                     
                     trt_type*beta_trt_baseline_OS+
                     
                     localization_0*beta_localization_OS +
                     localization_1*beta_localization_OS +
                     localization_2*beta_localization_OS + 
                     localization_3 *beta_localization_OS + 
                     localization_4*beta_localization_OS +
                     localization_5*beta_localization_OS +
                     
                     multifocal*beta_multifocal +
                     pTERT*beta_pTERT_OS +
                     tumor_size*beta_tumor_size_OS+
                     EGFR*beta_EGFR_OS +
                     g7l10*beta_g7l10_OS+
                     unmeas_binary*beta_unmeas_binary
  )
  
  lp_OS_z1 <- with(registry,
                   int_OS +
                     beta_M_OS*M + 
                     beta_B_OS*B +
                     beta_age_OS*(age - 60) +
                     beta_Unmeas_OS*unmeasured +
                     WHOm2*beta_WHOm2_OS +
                     sex*beta_sex_OS +
                     resection*beta_resection_OS +
                     trt_type*beta_trt_baseline_OS+
                     
                     localization_0*beta_localization_OS +
                     localization_1*beta_localization_OS +
                     localization_2*beta_localization_OS + 
                     localization_3 *beta_localization_OS + 
                     localization_4*beta_localization_OS +
                     localization_5*beta_localization_OS +                     
                     
                     multifocal*beta_multifocal +
                     pTERT*beta_pTERT_OS +
                     tumor_size*beta_tumor_size_OS+
                     EGFR*beta_EGFR_OS +
                     g7l10*beta_g7l10_OS +
                     unmeas_binary*beta_unmeas_binary+
                     
                     beta_treatment_OS 

                   
  )
  
  
  registry$lp_OS_z0<-lp_OS_z0
  registry$lp_OS_z1<-lp_OS_z1
  
  # time under control and time under treatment 
  
  registry$T_true_z0 <- rexp(N, rate = lambda0 * exp(lp_OS_z0))
  registry$T_true_z1 <- rexp(N, rate = lambda0 * exp(lp_OS_z1))
  
  
  registry$time_z0  <- pmin(registry$T_true_z0, C)
  registry$event_z0 <- as.numeric(registry$T_true_z0 <= C)
  
  registry$time_z1  <- pmin(registry$T_true_z1, C)
  registry$event_z1 <- as.numeric(registry$T_true_z1 <= C)
  
  registry |>
    dplyr::mutate(time = time_z0, event = event_z0,
                  resection=as.factor(resection),
                  trt_type=as.factor(trt_type)) |> 
    dplyr::select(-unmeas_binary) 
}
