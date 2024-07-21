library(tidyverse)
library(odin)
library(patchwork)
library(stringr)

# beta: Number of effective contacts made by a typical infectious individual per day
# mu: Probability of death due to infection 
# theta: 1 – vaccine efficacy against infection (VE_infection)/100%
# kappa: 1 – vaccine efficacy against death given infection (VE_(death|infection))/100%
# gamma: Recovery rate per day

# (1) Case 1: Time-invariant parameters ----------------------------------------
case1_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case1_time_invariant_par.R")

df_all_counterfactuals_case1 <-
  map_dfr(data.frame(
    "alpha==0.90" = 0.90,
    "alpha==0.70" = 0.70,
    "alpha==0.00" = 0
  ), function(par) {
    
    alpha <- par[1]

    model <- case1_generator$new(
      S_s_ini = 2e4 * (1 - alpha) - 20 * (1 - alpha),
      S_v_ini = 2e4 * alpha - 20 * alpha,
      I_s_ini = 20 * (1 - alpha),
      I_v_ini = 20 * alpha,
      R_s_ini = 0,
      R_v_ini = 0,
      D_s_ini = 0,
      D_v_ini = 0,
      Cum_vax_ini = 2e4 * alpha,
      Cum_inf_ini_v = 20 * alpha,
      Cum_inf_ini_s = 20 * (1 - alpha),
      beta = 0.15,  
      theta = 0.5,
      kappa = 0.1,
      gamma = 0.07,
      mu = 0.01
    )
    as_tibble(cbind(
      model$run(0:(730 * 1)),
      VE_infection = 1 - 0.5,
      alpha = alpha
    ))
  }, .id = "par")

# a function to collapse by vax status
collapse_vax_status <- function(dataframe){
    dataframe %>%
      mutate(
        S = S_s + S_v,
        I = I_s + I_v,
        D = D_s + D_v,
        R = R_s + R_v,
        N = S + I + R,
        N_v = S_v + I_v + R_v,
        N_s = S_s + I_s + R_s,
        Cum_inf = Cum_inf_s + Cum_inf_v
      )
}

df_all_counterfactuals_case1 <- collapse_vax_status(df_all_counterfactuals_case1)

# (2) Case 2: Increasing effective contacts (beta) -----------------------------
case2_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case2_inc_beta.R")
df_all_counterfactuals_case2 <-
  map_dfr(data.frame(
    "alpha==0.90" = 0.90,
    "alpha==0.70" = 0.70,
    "alpha==0.00" = 0
    ), function(par) {
      
      alpha <- par[1]

      flux_beta_t = c(0, 300)
      flux_beta_y = c(0.15, 0.6)
      
      model <- case2_generator$new(
        S_s_ini = 2e4 * (1-alpha) - 20 * (1-alpha),
        S_v_ini = 2e4 * alpha - 20 * alpha,
        I_s_ini = 20 * (1-alpha),
        I_v_ini = 20 * alpha,
        R_s_ini = 0,
        R_v_ini = 0,
        D_s_ini = 0,
        D_v_ini = 0,
        Cum_vax_ini = 2e4 * alpha,
        Cum_inf_ini_v = 20 * alpha,
        Cum_inf_ini_s = 20 * (1 - alpha),
        flux_beta_t = flux_beta_t,
        flux_beta_y = flux_beta_y,
        theta = 0.5,
        kappa = 0.1,
        gamma = 0.07,
        mu = 0.01) 
      
      as_tibble(cbind(model$run(0:(730 * 1)),
                      VE_infection = 1-0.5,
                      alpha = alpha))
    }, .id = "par")

df_all_counterfactuals_case2 <- collapse_vax_status(df_all_counterfactuals_case2)

# (3) Case 3: Increasing IFR (mu) ----------------------------------------------
case3_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case3_inc_mu.R")
df_all_counterfactuals_case3 <-
  map_dfr(data.frame(
    "alpha==0.90" = 0.90,
    "alpha==0.70" = 0.70,
    "alpha==0.00" = 0
  ), function(par) {
    
    alpha <- par[1]

    flux_mu_t = c(0, 300)
    flux_mu_y = c(0.01, 0.1)
    
    model <- case3_generator$new(
      S_s_ini = 2e4 * (1-alpha) - 20 * (1-alpha),
      S_v_ini = 2e4 * alpha - 20 * alpha,
      I_s_ini = 20 * (1-alpha),
      I_v_ini = 20 * alpha,
      R_s_ini = 0,
      R_v_ini = 0,
      D_s_ini = 0,
      D_v_ini = 0,
      Cum_vax_ini = 2e4 * alpha,
      Cum_inf_ini_v = 20 * alpha,
      Cum_inf_ini_s = 20 * (1 - alpha),
      beta = 0.15,
      flux_mu_t = flux_mu_t,
      flux_mu_y = flux_mu_y,
      theta = 0.5,
      kappa = 0.1,
      gamma = 0.07)
    
    as_tibble(cbind(model$run(0:(730 * 1)),
                    VE_infection = 1-0.5,
                    alpha = alpha))
    
  }, .id = "par")

df_all_counterfactuals_case3 <- collapse_vax_status(df_all_counterfactuals_case3)

# (4) Case 4: Waning VEs -------------------------------------------------------
case4_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case4_waning_VEs.R")
df_all_counterfactuals_case4 <-
  map_dfr(data.frame(
    "alpha==0.90" = 0.90,
    "alpha==0.70" = 0.70,
    "alpha==0.00" = 0
  ), function(par) {
    
    alpha <- par[1]

    flux_theta_waned_t = c(0, 100, 300, 730)
    flux_theta_waned_y = c(0.5, 0.5, 1, 1)
    
    flux_kappa_waned_t = c(0, 100, 300, 730)
    flux_kappa_waned_y = c(0.1, 0.1, 1, 1)
    
    model <- case4_generator$new(
      S_s_ini = 2e4 * (1-alpha) - 20 * (1-alpha),
      S_v_ini = 2e4 * alpha - 20 * alpha,
      I_s_ini = 20 * (1-alpha),
      I_v_ini = 20 * alpha,
      R_s_ini = 0,
      R_v_ini = 0,
      D_s_ini = 0,
      D_v_ini = 0,
      Cum_vax_ini = 2e4 * alpha,
      Cum_inf_ini_v = 20 * alpha,
      Cum_inf_ini_s = 20 * (1 - alpha),
      beta = 0.15,
      mu = 0.01,
      flux_theta_waned_t = flux_theta_waned_t,
      flux_theta_waned_y = flux_theta_waned_y, 
      flux_kappa_waned_t = flux_kappa_waned_t,
      flux_kappa_waned_y = flux_kappa_waned_y,
      gamma = 0.07)
    
    as_tibble(cbind(model$run(0:(730 * 1)),
                    VE_infection = 1-0.5,
                    alpha = alpha))
    
  }, .id = "par")

df_all_counterfactuals_case4 <- collapse_vax_status(df_all_counterfactuals_case4)

# (5) Case 5: Cases 2 and 4 ----------------------------------------------------
case5_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case5_inc_beta_waning_VEs.R")

df_all_counterfactuals_case5 <-
  map_dfr(data.frame(
    "alpha==0.90" = 0.90,
    "alpha==0.70" = 0.70,
    "alpha==0.00" = 0
  ), function(par) {
    
    alpha <- par[1]

    flux_theta_waned_t = c(0, 100, 300, 730)
    flux_theta_waned_y = c(0.5, 0.5, 1, 1)
    
    flux_kappa_waned_t = c(0, 100, 300, 730)
    flux_kappa_waned_y = c(0.1, 0.1, 1, 1)
    
    flux_beta_t = c(0, 300)
    flux_beta_y = c(0.15, 0.6)
    
    model <- case5_generator$new(
      S_s_ini = 2e4 * (1-alpha) - 20 * (1-alpha),
      S_v_ini = 2e4 * alpha - 20 * alpha,
      I_s_ini = 20 * (1-alpha),
      I_v_ini = 20 * alpha,
      R_s_ini = 0,
      R_v_ini = 0,
      D_s_ini = 0,
      D_v_ini = 0,
      Cum_vax_ini = 2e4 * alpha,
      Cum_inf_ini_v = 20 * alpha,
      Cum_inf_ini_s = 20 * (1 - alpha),
      flux_beta_t = flux_beta_t,
      flux_beta_y = flux_beta_y,
      mu = 0.01,
      flux_theta_waned_t = flux_theta_waned_t,
      flux_theta_waned_y = flux_theta_waned_y, 
      flux_kappa_waned_t = flux_kappa_waned_t,
      flux_kappa_waned_y = flux_kappa_waned_y,
      gamma = 0.07
      )
    
    as_tibble(cbind(model$run(0:(730 * 1)),
                    VE_infection = 1-0.5,
                    alpha = alpha))
    
  }, .id = "par")

df_all_counterfactuals_case5 <- collapse_vax_status(df_all_counterfactuals_case5)
