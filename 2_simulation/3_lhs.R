# Latin hypercube sampling 
library(tidyverse)
library(lhs)
library(ggforce)
set.seed(123)

# <I> Test whether Claim 1 holds under Scenario 1 ------------------------------
## (1) LHS ----------------------------------------------------------------------
ranges <- list(
  "R0" = c(1, 20),
  "alpha_observed" = c(0.01 ,1),
  "alpha_hypothetical" = c(0.01, 1),
  "mu" = c(0, 0.1),
  "gamma" = c(1 / 3, 1 / 14),
  "theta" = c(0, 1),
  "kappa" = c(0, 1),
  "T_VE_infection_unwaned" = c(0, 365),
  "T_VE_infection_waning" = c(0, 365),
  "T_VE_death_unwaned" = c(0, 365),
  "T_VE_death_waning" = c(0, 365)
)

ndraws <- 1000
parms_values <- as.data.frame(randomLHS(ndraws, 11)) 
names(parms_values) <- names(ranges)
for(parm in names(ranges)){
  parms_values[,parm] <- parms_values[,parm]*(ranges[[parm]][2]-ranges[[parm]][1]) + ranges[[parm]][1]
}

## (2) Test Claim 1 under Case 1 ------------------------------------------------
# call Odin model
case1_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case1_time_invariant_par.R")

# a function to run model under Scenario 1
simulate_case1 <- function(par, alpha_hypothetical, alpha_observed){

  # Simulate counterfactuals under (alpha_observed, alpha_hypothetical)
  out <- map_dfr(data.frame("alpha==alpha_observed" = alpha_observed,
                            "alpha==alpha_hypothetical" = alpha_hypothetical
                            ), function(prop_vax_var) {
                              
                              alpha <- prop_vax_var
                              
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
                                beta = as.numeric(par["R0"]) * as.numeric(par["gamma"]),
                                theta = as.numeric(par["theta"]),
                                kappa = as.numeric(par["kappa"]),
                                gamma = as.numeric(par["gamma"]),
                                mu = as.numeric(par["mu"])
                              )
                              as_tibble(cbind(
                                model$run(0:(730 * 1), method = "ode45"),
                                VE_infection = 1 - as.numeric(par["theta"]),
                                alpha = alpha,
                                alpha_observed = as.numeric(par["alpha_observed"]),
                                alpha_hypothetical = as.numeric(par["alpha_hypothetical"]),
                                abs_alpha_diff = abs(as.numeric(par["alpha_observed"]) - as.numeric(par["alpha_hypothetical"]))
                              ))
                            }, .id = "par")
  
  out <- collapse_vax_status(out)
  
  return(out)}

# a function to test Claim 1 under case 1
test_claim1_case1 <- function(par){
  
  # simulate
  out <- simulate_case1(par, alpha_hypothetical = par["alpha_hypothetical"], alpha_observed = par["alpha_observed"])
  
  # flag if the numerical integration is unsuccessful
  tt <- tryCatch(simulate_case1(par, alpha_hypothetical = par["alpha_hypothetical"], alpha_observed = par["alpha_observed"]),
                 error=function(e) e, 
                 warning=function(w) w)
  
  flag_unsuccessful_numerical_integration <- ifelse(is(tt,"warning"), "integration unsuccessful", 0)
  
  # get PDE(t, alpha_hypothetical, alpha_observed)
  PDE <- get_PDE(out, alpha_hypothetical = par["alpha_hypothetical"], alpha_observed = par["alpha_observed"])
  
  # get POE
  POE <- get_POE(out)

  # Evaluate
  POE_PDE <- PDE %>%
    left_join(POE, by = "t") 
  
  days_POE_less_PDE_infection <- sum(POE_PDE$POE_infection - POE_PDE$PDE_infection < -10e-7)
  days_POE_less_PDE_death <- sum(POE_PDE$POE_death - POE_PDE$PDE_death < -10e-7)
  
  return(c(days_POE_less_PDE_infection=days_POE_less_PDE_infection,
           days_POE_less_PDE_death=days_POE_less_PDE_death,
           flag_unsuccessful_numerical_integration))
}

## (2.1) Total days in which Claim 1 does not hold under Scenario 1
total_days_claim1_is_false_case1 <- apply(parms_values, 1, test_claim1_case1)
colnames(total_days_claim1_is_false_case1) <- 1:ndraws

## (2.2) Check no par gives unsuccessful numerical integration
(no_integral_par_case1 <- which(total_days_claim1_is_false_case1[3,]!=0, TRUE))

## (2.3) Check none of the LHS samples gives POE < PDE
(claim1a_is_false_case1_infection <- which(total_days_claim1_is_false_case1[1,] > 0, TRUE))
(claim1a_is_false_case1_death <- which(total_days_claim1_is_false_case1[2,] > 0, TRUE))

# <II> Test whether Claim 1a holds under Scenario 4 ----------------------------
## (1) Test Claim 1a under Scenario 4 --------------------------------------------
# call Odin model
case4_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case4_waning_VEs.R")

# a function to run model under Scenario 4
simulate_case4 <- function(par, alpha_hypothetical=0, alpha_observed){
  
  flux_theta_waned_t <- c(0,
                          round(par["T_VE_infection_unwaned"]),
                          round(par["T_VE_infection_unwaned"]) + round(par["T_VE_infection_waning"]),
                          730)
  flux_theta_waned_y <- c(par["theta"], par["theta"], 1, 1)
  
  flux_kappa_waned_t <- c(0,
                          round(par["T_VE_death_unwaned"]),
                          round(par["T_VE_death_unwaned"]) + round(par["T_VE_death_waning"]),
                          730)
  flux_kappa_waned_y <- c(par["kappa"], par["kappa"], 1, 1)
  
  # Simulate counterfactuals under (alpha_hypothetical = 0, alpha_observed)
  out <- map_dfr(data.frame("alpha==alpha_observed" = alpha_observed,
                            "alpha==0" = 0), function(prop_vax_var) {
                              
                              alpha <- prop_vax_var
                              
                              # Assign model par
                              model <- case4_generator$new(
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
                                beta = as.numeric(par["R0"]) * as.numeric(par["gamma"]),
                                flux_theta_waned_t = as.numeric(flux_theta_waned_t),
                                flux_theta_waned_y = as.numeric(flux_theta_waned_y),
                                flux_kappa_waned_t = as.numeric(flux_kappa_waned_t),
                                flux_kappa_waned_y = as.numeric(flux_kappa_waned_y),
                                gamma = as.numeric(par["gamma"]),
                                mu = as.numeric(par["mu"])
                              )
                              
                              # Run model
                              as_tibble(cbind(model$run(0:(730 * 1), method = "ode45"),
                                              VE_infection = 1 - as.numeric(par["theta"]),
                                              alpha = alpha,
                                              alpha_observed = as.numeric(par["alpha_observed"]),
                                              alpha_hypothetical = 0,
                                              abs_alpha_diff = abs(as.numeric(par["alpha_observed"]) - 0)))
                            }, .id = "par")
  
  out <- collapse_vax_status(out)
  
  return(out)
}

# a function to test Claim 2 under case 4
test_claim1a_case4 <- function(par){
  
  # simulate
  out <- simulate_case4(par, alpha_hypothetical = 0, alpha_observed = par["alpha_observed"])
  
  # flag if the numerical integration is unsuccessful
  tt <- tryCatch(simulate_case4(par, alpha_hypothetical = par["alpha_hypothetical"], alpha_observed = par["alpha_observed"]),
                 error=function(e) e, 
                 warning=function(w) w)
  
  flag_unsuccessful_numerical_integration <- ifelse(is(tt,"warning"), "integration unsuccessful", 0)
  
  # get PDE(t, alpha_hypothetical=0, alpha_observed)
  PDE <- get_PDE(out, alpha_hypothetical = 0, alpha_observed = as.numeric(par["alpha_observed"]))
  
  # get POE
  POE <- get_POE(out)
  
  # Evaluate
  POE_PDE <- PDE %>%
    left_join(POE, by = "t") 
  
  days_POE_less_PDE_infection <- sum(POE_PDE$POE_infection - POE_PDE$PDE_infection < -10e-7)
  days_POE_less_PDE_death <- sum(POE_PDE$POE_death - POE_PDE$PDE_death < -10e-7)
  
  return(c(days_POE_less_PDE_infection = days_POE_less_PDE_infection,
           days_POE_less_PDE_death = days_POE_less_PDE_death,
           flag = flag_unsuccessful_numerical_integration))
  
}

## (2.1) Total days in which Claim 1a does not hold under Scenario 4
total_days_claim1a_is_false_case4 <- apply(parms_values, 1, test_claim1a_case4)
colnames(total_days_claim1a_is_false_case4) <- 1:ndraws

## (2.2) Check no par gives unsuccessful numerical integration
(no_integral_par_case4 <- which(total_days_claim1a_is_false_case4[3,]!=0, TRUE))

## (2.3) Check none of the LHS samples gives POE < PDE
(claim1a_is_false_case4_infection <- which(total_days_claim1a_is_false_case4[1,] > 0, TRUE))
(claim1a_is_false_case4_death <- which(total_days_claim1a_is_false_case4[2,] > 0, TRUE))

## <III> Randomly draw 50 samples from LHS and plot PDEs and POEs ------------------
### (3.1) randomly draw samples
set.seed(123)
par_id <- 1:ndraws 
subsamples_id <- sample(par_id, size=50, replace = FALSE)
subsamples <- parms_values[subsamples_id, ]
subsamples <- cbind(subsamples, par_id=subsamples_id)

### (3.2) get PDEs and POEs for plotting
# a function to get PDE and POE for plotting under Case 1
get_POE_PDE_for_plotting <- function(subpar){
  
  # simulate case 1 under claim 1
  out_case1 <- simulate_case1(subpar, alpha_hypothetical = as.numeric(subpar["alpha_hypothetical"]), alpha_observed = as.numeric(subpar["alpha_observed"]))
  
  # get PDE(t, alpha_hypothetical, alpha_observed)
  PDE_case1 <- get_PDE(out_case1, alpha_hypothetical = as.numeric(subpar["alpha_hypothetical"]), alpha_observed = as.numeric(subpar["alpha_observed"]))
  
  # get POE
  POE_case1 <- get_POE(out_case1)
  
  # Left join
  POE_PDE_case1 <- PDE_case1 %>%
    left_join(POE_case1, by = "t") %>%
    mutate(case = "Scenario 1") %>%
    select(par, t, alpha, alpha_hypothetical, alpha_observed, PDE_infection, PDE_death, POE_infection, POE_death, case)
  
  # simulate case 4 under claim 1a
  out_case4 <- simulate_case4(subpar, alpha_hypothetical = subpar["alpha_hypothetical"], alpha_observed = subpar["alpha_observed"])
  
  # get PDE(t, alpha_hypothetical, alpha_observed)
  PDE_case4 <- get_PDE(out_case4, alpha_hypothetical = 0, alpha_observed = as.numeric(subpar["alpha_observed"]))
  
  # get POE
  POE_case4 <- get_POE(out_case4)
  
  # Left join
  POE_PDE_case4 <- PDE_case4 %>%
    left_join(POE_case4, by = "t") %>%
    mutate(case = "Scenario 4") %>%
    select(par, t, alpha, alpha_hypothetical, alpha_observed, PDE_infection, PDE_death, POE_infection, POE_death, case)
  
  POE_PDE <- rbind(POE_PDE_case1, POE_PDE_case4)
  
  return(POE_PDE)
  
}

### (3.3) plot
df_POE_PDE_for_plotting <- apply(subsamples, 1, get_POE_PDE_for_plotting) 

df_POE_PDE_for_plotting <- do.call(rbind, df_POE_PDE_for_plotting)  

df_POE_PDE_for_plotting_long <- df_POE_PDE_for_plotting %>% 
  select(case, alpha, t, PDE_infection, PDE_death, POE_infection, POE_death) %>% 
  pivot_longer(
    cols = -c("case","alpha", "t"),
    names_to = c("estimand", "outcome"),
    names_sep = "_",
    values_to = "count"
  )

subsamples_POE_PDE_death <- ggplot() +
  geom_line(data = df_POE_PDE_for_plotting_long %>% filter(outcome=="death"), aes(x = t, y= count, group = alpha), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y= bquote(atop('POE'^'death'*~or
                      ~'PDE'^'death')),
       x="Day") +
  facet_grid(case ~ estimand, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

subsamples_POE_PDE_infection <- ggplot() +
  geom_line(data = df_POE_PDE_for_plotting_long %>% filter(outcome=="infection"), aes(x = t, y= count, group = alpha), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y= bquote(atop('POE'^'infection'*~or
                      ~'PDE'^'infection')),
       x="Day") +
  facet_grid(case ~ estimand, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

# collate the plots
col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_eFig3 <- "
fffffffggggggg
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii
hhhhhhhiiiiiii"

plotlist_eFig3 <-
  list(
    f = col1,
    g = col2,
    h = subsamples_POE_PDE_infection,
    i = subsamples_POE_PDE_death)

plot_subsamples_POE_PDE_labeled <- wrap_plots(plotlist_eFig3, guides = 'collect', design = layoutplot_eFig3) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig3.png", plot_subsamples_POE_PDE_labeled, width = 11, height=6, dpi=300, units="in")
