# Latin hypercube sampling 
library(lhs)
set.seed(123)

# <I> LHS ----------------------------------------------------------------------
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

# <II> functions ---------------------------------------------------------------
## (1) a function to simulate the claims ---------------------------------------
simulate_claims <-
  function(par, claim_scenario) {
    
    par <- unlist(par)
    
    if(claim_scenario=="claim1_scenario1"){
      
      # Claim 1 applies to any alpha, alpha'
      ALPHA_OBSERVED <- par["alpha_observed"]
      ALPHA_HYPOTHETICAL <- par["alpha_hypothetical"]
      
      # Scenario 1 has time-invariant parameters
      flux_theta_y <- rep(par["theta"], 4)
      flux_kappa_y <- rep(par["kappa"], 4)
      
    }else{
      # Claim 1a applies to alpha=0, alpha' in [0,1] 
      ALPHA_OBSERVED <- par["alpha_observed"]
      ALPHA_HYPOTHETICAL <- 0
      
      # Scenario 4 has time-invariant parameters
      flux_theta_y <- c(par["theta"], par["theta"], 1, 1)
      flux_kappa_y <- c(par["kappa"], par["kappa"], 1, 1)
    }
    
    out_scenario <-
      map_dfr(data.frame(
        "alpha==alpha_observed" = ALPHA_OBSERVED,
        "alpha==alpha_hypothetical" = ALPHA_HYPOTHETICAL), function(prop_vax_var) {
        
        alpha <- prop_vax_var
        
        model <- model_generator$new(
          S_s_ini = (2e4 - 20) * (1 - alpha),
          S_v_ini = (2e4 - 20) * alpha,
          I_s_ini = 20 * (1 - alpha),
          I_v_ini = 20 * alpha,
          R_s_ini = 0,
          R_v_ini = 0,
          D_s_ini = 0,
          D_v_ini = 0,
          Cum_vax_ini = 2e4 * alpha,
          Cum_inf_ini_v = 20 * alpha,
          Cum_inf_ini_s = 20 * (1 - alpha),
          flux_beta_t = c(0, 300),
          # no increasing beta
          flux_beta_y = rep(par["R0"] * par["gamma"], 2),
          flux_mu_t = c(0, 300),
          # no increasing mu
          flux_mu_y = rep(par["mu"], 2),
          flux_theta_t = c(0,
                           round(par["T_VE_infection_unwaned"]),
                           round(par["T_VE_infection_unwaned"]) + round(par["T_VE_infection_waning"]),
                           730),
          # no waning VE infection
          flux_theta_y = flux_theta_y,
          flux_kappa_t = c(0,
                           round(par["T_VE_death_unwaned"]),
                           round(par["T_VE_death_unwaned"]) + round(par["T_VE_death_waning"]),
                           730),
          # no waning VE death
          flux_kappa_y = flux_kappa_y,
          gamma = par["gamma"]
        )
        
        out_scenario <- as_tibble(
          cbind(
            model$run(0:(730 * 1), method = "ode45"),
            alpha = alpha,
            alpha_observed = ALPHA_OBSERVED,
            alpha_hypothetical = ALPHA_HYPOTHETICAL)
        )
        
        out_scenario <-
          collapse_vax_status(out_scenario)
      }, .id = "par")
    
    return(out_scenario)
  }

## (2) a function to test claims -----------------------------------------------
test_claims <- function(par, claim_scenario){
  
  par <- unlist(par)
  
  if(claim_scenario=="claim1_scenario1"){
    # Claim 1 applies to any alpha, alpha'
    ALPHA_OBSERVED <- par["alpha_observed"]
    ALPHA_HYPOTHETICAL <- par["alpha_hypothetical"]
    
  }else{
    # Claim 1a applies to alpha=0, alpha' in [0,1] 
    ALPHA_OBSERVED <- par["alpha_observed"]
    ALPHA_HYPOTHETICAL <- 0
    
  }
  
  # simulate
  out_scenario <- simulate_claims(par, claim_scenario)
  
  # flag if the numerical integration is unsuccessful
  tt <- tryCatch(simulate_claims(par, claim_scenario),
                 error=function(e) e, 
                 warning=function(w) w)
  
  flag_unsuccessful_numerical_integration <- ifelse(is(tt,"warning"), "integration unsuccessful", 0)
  
  # evaluate
  PDE_POE <- get_PDE_POE(out_scenario, alpha_hypothetical = ALPHA_HYPOTHETICAL, alpha_observed = ALPHA_OBSERVED)
  
  days_POE_less_PDE_infection <- sum(PDE_POE$POE_infection - PDE_POE$PDE_infection < -10e-7)
  days_POE_less_PDE_death <- sum(PDE_POE$POE_death - PDE_POE$PDE_death < -10e-7)
  
  return(c(days_POE_less_PDE_infection=days_POE_less_PDE_infection,
           days_POE_less_PDE_death=days_POE_less_PDE_death,
           flag=flag_unsuccessful_numerical_integration))
}

# <III> test whether Claim 1 holds under Scenario 1 ----------------------------
## (1) total days in which Claim 1 does not hold under Scenario 1
total_days_claim1_is_false_scenario_1 <- apply(parms_values, 1, function(par_vec) test_claims(par_vec, claim_scenario = "claim1_scenario1"))
colnames(total_days_claim1_is_false_scenario_1) <- 1:ndraws

## (2) check no par gives unsuccessful numerical integration
(no_integral_par_scenario_1 <- which(total_days_claim1_is_false_scenario_1["flag",]!=0, TRUE))

## (3) check none of the LHS samples gives POE < PDE
(claim1a_is_false_scenario_1_infection <- which(total_days_claim1_is_false_scenario_1["days_POE_less_PDE_infection",] > 0, TRUE))
(claim1a_is_false_scenario_1_death <- which(total_days_claim1_is_false_scenario_1["days_POE_less_PDE_death",] > 0, TRUE))

# <IV> test whether Claim 1a holds under Scenario 4 ----------------------------
## (1) total days in which Claim 1a does not hold under Scenario 4
total_days_claim1a_is_false_scenario_4 <- apply(parms_values, 1, function(par_vec) test_claims(par_vec, claim_scenario = "claim1a_scenario4"))
colnames(total_days_claim1a_is_false_scenario_4) <- 1:ndraws

## (2) Check no par gives unsuccessful numerical integration
(no_integral_par_scenario_4 <- which(total_days_claim1a_is_false_scenario_4["flag",]!=0, TRUE))

## (3) Check none of the LHS samples gives POE < PDE
(claim1a_is_false_scenario_4_infection <- which(total_days_claim1a_is_false_scenario_4["days_POE_less_PDE_infection",] > 0, TRUE))
(claim1a_is_false_scenario_4_death <- which(total_days_claim1a_is_false_scenario_4["days_POE_less_PDE_death",] > 0, TRUE))

# <V> Randomly draw 50 samples from LHS and plot PDEs and POEs ---------------
## (1) randomly draw samples -------------------------------------------------
set.seed(123)
par_id <- 1:ndraws 
subsamples_id <- sample(par_id, size=50, replace = FALSE)
subsamples <- parms_values[subsamples_id, ]
subsamples <- cbind(subsamples, par_id=subsamples_id)

## (2) get PDEs and POEs for plotting ----------------------------------------
# a function to get PDE and POE for plotting 
get_PDE_POE_for_plotting <- function(subpar){
  
  subpar <- unlist(subpar)
  
  ALPHA_OBSERVED <- subpar["alpha_observed"]
  ALPHA_HYPOTHETICAL <- subpar["alpha_hypothetical"]
  
  # simulate Scenario 1 under claim 1
  out_scenario_1 <- simulate_claims(subpar, claim_scenario="claim1_scenario1")
  
  # get PDE and POE
  PDE_POE_scenario_1 <- get_PDE_POE(out_scenario_1, alpha_hypothetical = ALPHA_HYPOTHETICAL, alpha_observed = ALPHA_OBSERVED) %>%
    mutate(scenario = "Scenario 1",
           subpar_id = subpar["par_id"]) 
  
  # simulate scenario 4 under claim 1a
  out_scenario_4 <- simulate_claims(subpar, claim_scenario="claim1a_scenario4")
  
  # get PDE and POE
  PDE_POE_scenario_4 <- get_PDE_POE(out_scenario_4, alpha_hypothetical = 0, alpha_observed = ALPHA_OBSERVED) %>%
    mutate(scenario = "Scenario 4",
           subpar_id = subpar["par_id"]) 
  
  PDE_POE <- rbind(PDE_POE_scenario_1, PDE_POE_scenario_4)
  return(PDE_POE)
}

### (3) plot
df_PDE_POE_for_plotting <- apply(subsamples, 1, get_PDE_POE_for_plotting) 

df_PDE_POE_for_plotting <- do.call(rbind, df_PDE_POE_for_plotting)  

df_PDE_POE_for_plotting_long <- df_PDE_POE_for_plotting %>% 
  select(scenario, subpar_id, t, PDE_infection, PDE_death, POE_infection, POE_death) %>% 
  pivot_longer(
    cols = -c("scenario","subpar_id", "t"),
    names_to = c("estimand", "outcome"),
    names_sep = "_",
    values_to = "count"
  )

subsamples_PDE_POE_death <- ggplot() +
  geom_line(data = df_PDE_POE_for_plotting_long %>% filter(outcome=="death"), aes(x = t, y= count, group = subpar_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y= bquote(atop('POE'^'death'*~or
                      ~'PDE'^'death')),
       x="Day") +
  facet_grid(scenario ~ estimand, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

subsamples_PDE_POE_infection <- ggplot() +
  geom_line(data = df_PDE_POE_for_plotting_long %>% filter(outcome=="infection"), aes(x = t, y= count, group = subpar_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y= bquote(atop('POE'^'infection'*~or
                      ~'PDE'^'infection')),
       x="Day") +
  facet_grid(scenario ~ estimand, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

# collate the plots
col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_eFig2 <- "
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

plotlist_eFig2 <-
  list(
    f = col1,
    g = col2,
    h = subsamples_PDE_POE_infection,
    i = subsamples_PDE_POE_death)

plot_subsamples_PDE_POE_labeled <- wrap_plots(plotlist_eFig2, guides = 'collect', design = layoutplot_eFig2) 

ggsave("3_figures/eFig2.png", plot_subsamples_PDE_POE_labeled, width = 11, height=6, dpi=300, units="in")
