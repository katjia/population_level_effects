library(tidyverse)
library(odin)
library(patchwork)
library(stringr)

# beta: The mumber of effective contacts made by a typical infectious individual per day
# mu: Probability of death due to infection 
# theta: 1 – vaccine efficacy against infection (VE_infection)/100%
# kappa: 1 – vaccine efficacy against death given infection (VE_(death|infection))/100%
# gamma: Recovery rate per day

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

# <I> compile Odin model -------------------------------------------------------
model_generator <- odin::odin("1_model/Scenarios_1_to_5.R")

# <II> assign model parameters (eTable 1) --------------------------------------
  ## Assign time-invariant values as default
  flux_beta_t <- matrix(rep(c(0, 300), 5), nrow=5, byrow=TRUE)
  flux_beta_y <- matrix(rep(c(0.15, 0.15), 5), nrow=5, byrow=TRUE)
  
  flux_mu_t <- matrix(rep(c(0, 300), 5), nrow=5, byrow=TRUE)
  flux_mu_y <- matrix(rep(c(0.01, 0.01), 5), nrow=5, byrow=TRUE)
  
  flux_theta_t <- matrix(rep(c(0, 100, 300, 730), 5), nrow=5, byrow=TRUE) 
  flux_theta_y <- matrix(rep(c(0.5, 0.5, 0.5, 0.5), 5), nrow=5, byrow=TRUE) 
  
  flux_kappa_t <- matrix(rep(c(0, 100, 300, 730), 5), nrow=5, byrow=TRUE) 
  flux_kappa_y <- matrix(rep(c(0.1, 0.1, 0.1, 0.1), 5), nrow=5, byrow=TRUE) 
  
## (1) Scenario 1: Time-invariant parameters ------------------------------------
  ## The default (e.g., flux_beta_y[1, ] is time-invariant beta)

## (2) Scenario 2: Increasing effective contacts (beta) -------------------------
  flux_beta_t[2, ] = c(0, 300)
  flux_beta_y[2, ] = c(0.15, 0.6)

## (3) Scenario 3: Increasing IFR (mu) -------------------------------------
  flux_mu_t[3, ] = c(0, 300)
  flux_mu_y[3, ] = c(0.01, 0.1)

## (4) Scenario 4: Waning VEs ----------------------------------------------
  flux_theta_t[4, ] = c(0, 100, 300, 730)
  flux_theta_y[4, ] = c(0.5, 0.5, 1, 1)
  
  flux_kappa_t[4, ] = c(0, 100, 300, 730)
  flux_kappa_y[4, ] = c(0.1, 0.1, 1, 1)  

## (5) Scenario 5: Combination of Scenarios 2 and 4 ------------------------
  flux_beta_t[5,] = c(0, 300)
  flux_beta_y[5,] = c(0.15, 0.6)
  
  flux_theta_t[5,] = c(0, 100, 300, 730)
  flux_theta_y[5,] = c(0.5, 0.5, 1, 1)
  
  flux_kappa_t[5,] = c(0, 100, 300, 730)
  flux_kappa_y[5,] = c(0.1, 0.1, 1, 1)

# <III> run models under different alpha (proportion vaccinated) ---------------
  
  # A list for all 5 Scenarios
  df_all_counterfactuals_scenarios_1_to_5 <- list()
  
  for(scenario in 1:5){

  df_all_counterfactuals_scenario <-
    map_dfr(data.frame(
      "alpha==0.90" = 0.90,
      "alpha==0.70" = 0.70,
      "alpha==0.00" = 0
    ), function(par) {
      
      alpha <- par[1]
      
      model <- model_generator$new(
        S_s_ini = (2e4 - 20) * (1-alpha),
        S_v_ini = (2e4 - 20) * alpha,
        I_s_ini = 20 * (1-alpha),
        I_v_ini = 20 * alpha,
        R_s_ini = 0,
        R_v_ini = 0,
        D_s_ini = 0,
        D_v_ini = 0,
        Cum_vax_ini = 2e4 * alpha,
        Cum_inf_ini_v = 20 * alpha,
        Cum_inf_ini_s = 20 * (1 - alpha),
        flux_beta_t = flux_beta_t[scenario, ],
        flux_beta_y = flux_beta_y[scenario, ],
        flux_mu_t = flux_mu_t[scenario, ],
        flux_mu_y = flux_mu_y[scenario, ],
        flux_theta_t = flux_theta_t[scenario, ],
        flux_theta_y = flux_theta_y[scenario, ], 
        flux_kappa_t = flux_kappa_t[scenario, ],
        flux_kappa_y = flux_kappa_y[scenario, ],
        gamma = 0.07
      )
      
      out <- as_tibble(cbind(model$run(0:(730 * 1)),
                             alpha = alpha))
      
      collapse_vax_status(out)
      
    }, .id = "par")
  
  # compile a list
  name <- paste0("scenario_", scenario)
  df_all_counterfactuals_scenarios_1_to_5[[name]] <- df_all_counterfactuals_scenario

}
