## LHS to test whether Claim 1 holds under waning VEs
library(lhs)
library(ggforce)
set.seed(123)

# (1) LHS ----------------------------------------------------------------------
variable.parmranges <- list("R0"=c(1, 20), 
                            "alpha"=c(0.01, 1), 
                            "mu"=c(0, 0.1),
                            "gamma"=c(1/3, 1/14),
                            "theta"=c(0, 1), 
                            "kappa"=c(0, 1), 
                            "T_VE_infection_unwaned"=c(0, 365), 
                            "T_VE_infection_waning"=c(0, 365), 
                            "T_VE_death_unwaned"=c(0, 365), 
                            "T_VE_death_waning"=c(0, 365))

ndraws <- 1000
parmsdf.var <- as.data.frame(randomLHS(ndraws, 10)) 
names(parmsdf.var) <- names(variable.parmranges)
for(parm in names(variable.parmranges)){
  parmsdf.var[,parm] <- parmsdf.var[,parm]*(variable.parmranges[[parm]][2]-variable.parmranges[[parm]][1]) + variable.parmranges[[parm]][1]
}

# (2) Simulating traj to test Claim 1 under Case 1 (time-invariant parameters) ------------
test_claim1_case1 <- function(par){
  
  # Take par values from LHS
  prop_vax <- as.numeric(par["alpha"])
  
  # Simualte counterfactuals
  out <- map_dfr(data.frame("alpha=alpha1" = prop_vax,
                            "alpha=alpha0" = 0), function(prop_vax_var) {
                              
                              alpha <- prop_vax_var
                              
                              # Assign model par
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
                                theta = as.numeric(par["theta"]),
                                kappa = as.numeric(par["kappa"]),
                                beta = as.numeric(par["R0"]) * as.numeric(par["gamma"]),
                                gamma = as.numeric(par["gamma"]),
                                mu = as.numeric(par["mu"])
                              )
                              
                              # Run model
                              as_tibble(cbind(model$run(0:(730 * 1)),
                                              alpha = alpha))
                            }, .id = "par")
  
  ## Compute pY and pD among susc
  pY_pD <- out %>%
    select(alpha, t, Cum_inf_s, D_s) %>%
    mutate(
      pY_s = case_when(alpha > 0 ~ Cum_inf_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
      pD_s = case_when(alpha > 0 ~ D_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ D_s / 2e4)
    )
  
  ## Compute IE infection 
  IE_infection <- pY_pD %>%
    select(alpha, t, pY_s) %>%
    group_by(t) %>%
    summarise(IE_infection = last(pY_s) - first(pY_s))
  
  days_neg_IE_infection <- sum(IE_infection$IE_infection < -0.5/2e4)
  
  ## Compute IE death
  IE_death <- pY_pD %>%
    select(alpha, t, pD_s) %>%
    group_by(t) %>%
    summarise(IE_death = last(pD_s) - first(pD_s))
  
  days_neg_IE_death <- sum(IE_death$IE_death < -0.5/2e4)
  
  return(c(days_neg_IE_infection=days_neg_IE_infection,
           days_neg_IE_death=days_neg_IE_death))
}

## (2.1) Total days with neg IE
case1_days_neg_IE <- apply(parmsdf.var, 1, test_claim1_case4)
colnames(case1_days_neg_IE) <- 1:ndraws

## (2.2) Identify par combo for which IE is negative
(case1_neg_IE_infection_par <- which(case1_days_neg_IE[1,] > 0, TRUE))
(case1_neg_IE_death_par <- which(case1_days_neg_IE[2,] > 0, TRUE))

# (3) Simulating traj to test Claim 1 under Case 4 (waning VEs and time-invariant parameters) --------------
test_claim1_case4 <- function(par){

  # Take par values from LHS
  prop_vax <- as.numeric(par["alpha"])
  flux_theta_waned_t <-
    c(0,
      round(par["T_VE_infection_unwaned"]),
      round(par["T_VE_infection_unwaned"]) + round(par["T_VE_infection_waning"]),
      730)
  flux_theta_waned_y <- c(par["theta"], par["theta"], 1, 1)
  
  flux_kappa_waned_t <- 
    c(0,
      round(par["T_VE_death_unwaned"]),
      round(par["T_VE_death_unwaned"]) + round(par["T_VE_death_waning"]),
      730)
  flux_kappa_waned_y = c(par["kappa"], par["kappa"], 1, 1)
  
  # Simualte counterfactuals
  out <- map_dfr(data.frame("alpha=alpha1" = prop_vax,
                            "alpha=alpha0" = 0), function(prop_vax_var) {
                              
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
                                flux_theta_waned_t = as.numeric(flux_theta_waned_t),
                                flux_theta_waned_y = as.numeric(flux_theta_waned_y),
                                flux_kappa_waned_t = as.numeric(flux_kappa_waned_t),
                                flux_kappa_waned_y = as.numeric(flux_kappa_waned_y),
                                beta = as.numeric(par["R0"]) * as.numeric(par["gamma"]),
                                gamma = as.numeric(par["gamma"]),
                                mu = as.numeric(par["mu"])
                              )
                              
                              # Run model
                              as_tibble(cbind(model$run(0:(730 * 1)),
                                              alpha = alpha))
                            }, .id = "par")
      
  ## Compute pY and pD among susc
  pY_pD <- out %>%
    select(alpha, t, Cum_inf_s, D_s) %>%
    mutate(
      pY_s = case_when(alpha > 0 ~ Cum_inf_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
      pD_s = case_when(alpha > 0 ~ D_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ D_s / 2e4)
    )
      
  ## Compute IE infection 
  IE_infection <- pY_pD %>%
    select(alpha, t, pY_s) %>%
    group_by(t) %>%
    summarise(IE_infection = last(pY_s) - first(pY_s))
  
  days_neg_IE_infection <- sum(IE_infection$IE_infection < -0.5/2e4)
  
  ## Compute IE death
  IE_death <- pY_pD %>%
    select(alpha, t, pD_s) %>%
    group_by(t) %>%
    summarise(IE_death = last(pD_s) - first(pD_s))
  
  days_neg_IE_death <- sum(IE_death$IE_death < -0.5/2e4)
  
  return(c(days_neg_IE_infection=days_neg_IE_infection,
           days_neg_IE_death=days_neg_IE_death))
}

## (3.1) Total days with neg IE
case4_days_neg_IE <- apply(parmsdf.var, 1, test_claim1_case4)
colnames(case4_days_neg_IE) <- 1:ndraws

## (3.2) Identify par combo for which IE is negative
(case4_neg_IE_infection_par <- which(case4_days_neg_IE[1,] > 0, TRUE))
(case4_neg_IE_death_par <- which(case4_days_neg_IE[2,] > 0, TRUE))

# (4) Randomly draw n samples from LHS and plot the IEs ------------------------
## (4.1) randomly draw samples
set.seed(123)
subsamples_id <- sample(1:ndraws, size=50, replace = FALSE)
subsamples <- parmsdf.var[subsamples_id, ]
subsamples <- cbind(subsamples, par_id=subsamples_id)

## (4.2) get IEs for ploting
### Case 1
case1_get_IEs_for_plotting <- function(par) {
  
  # Take par values from LHS
  prop_vax <- as.numeric(par["alpha"])
  
  # Simualte counterfactuals
  out <- map_dfr(data.frame("alpha=alpha1" = prop_vax,
                            "alpha=alpha0" = 0), function(prop_vax_var) {
                              alpha <- prop_vax_var
                              
                              # Assign model par
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
                                theta = as.numeric(par["theta"]),
                                kappa = as.numeric(par["kappa"]),
                                beta = as.numeric(par["R0"]) * as.numeric(par["gamma"]),
                                gamma = as.numeric(par["gamma"]),
                                mu = as.numeric(par["mu"])
                              )
                              
                              # Run model
                              as_tibble(cbind(model$run(0:(730 * 1)),
                                              alpha = alpha))
                            }, .id = "par")
  
  ## Compute pY and pD among susc
  pY_pD <- out %>%
    select(alpha, t, Cum_inf_s, D_s) %>%
    mutate(
      pY_s = case_when(alpha > 0 ~ Cum_inf_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
      pD_s = case_when(alpha > 0 ~ D_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ D_s / 2e4)
    )
  
  ## Compute IE infection
  IE_infection <- pY_pD %>%
    select(alpha, t, pY_s) %>%
    group_by(t) %>%
    summarise(IE_infection = last(pY_s) - first(pY_s))
  
  ## Compute IE death
  IE_death <- pY_pD %>%
    select(alpha, t, pD_s) %>%
    group_by(t) %>%
    summarise(IE_death = last(pD_s) - first(pD_s))
  
  IEs <- IE_death %>%
    left_join(IE_infection, by="t") %>%
    cbind(alpha1=prop_vax,
          par_id=par["par_id"],
          row.names = NULL)
  
  return(IEs)
}

### Case 4
case4_get_IEs_for_plotting <- function(par) {
  
  # Take par values from LHS
  prop_vax <- as.numeric(par["alpha"])
  flux_theta_waned_t <-
    c(0,
      round(par["T_VE_infection_unwaned"]),
      round(par["T_VE_infection_unwaned"]) + round(par["T_VE_infection_waning"]),
      730)
  flux_theta_waned_y <- c(par["theta"], par["theta"], 1, 1)
  
  flux_kappa_waned_t <-
    c(0,
      round(par["T_VE_death_unwaned"]),
      round(par["T_VE_death_unwaned"]) + round(par["T_VE_death_waning"]),
      730)
  flux_kappa_waned_y = c(par["kappa"], par["kappa"], 1, 1)
  
  # Simualte counterfactuals
  out <- map_dfr(data.frame("alpha=alpha1" = prop_vax,
                            "alpha=alpha0" = 0), function(prop_vax_var) {
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
                                flux_theta_waned_t = as.numeric(flux_theta_waned_t),
                                flux_theta_waned_y = as.numeric(flux_theta_waned_y),
                                flux_kappa_waned_t = as.numeric(flux_kappa_waned_t),
                                flux_kappa_waned_y = as.numeric(flux_kappa_waned_y),
                                beta = as.numeric(par["R0"]) * as.numeric(par["gamma"]),
                                gamma = as.numeric(par["gamma"]),
                                mu = as.numeric(par["mu"])
                              )
                              
                              # Run model
                              as_tibble(cbind(model$run(0:(730 * 1)),
                                              alpha = alpha))
                            }, .id = "par")
  
  ## Compute pY and pD among susc
  pY_pD <- out %>%
    select(alpha, t, Cum_inf_s, D_s) %>%
    mutate(
      pY_s = case_when(alpha > 0 ~ Cum_inf_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
      pD_s = case_when(alpha > 0 ~ D_s / (2e4 * (1 - alpha)),
                       alpha == 0.0 ~ D_s / 2e4)
    )
  
  ## Compute IE infection
  IE_infection <- pY_pD %>%
    select(alpha, t, pY_s) %>%
    group_by(t) %>%
    summarise(IE_infection = last(pY_s) - first(pY_s))
  
  ## Compute IE death
  IE_death <- pY_pD %>%
    select(alpha, t, pD_s) %>%
    group_by(t) %>%
    summarise(IE_death = last(pD_s) - first(pD_s))
  
  IEs <- IE_death %>%
    left_join(IE_infection, by="t") %>%
    cbind(alpha1=prop_vax,
          par_id=par["par_id"],
          row.names = NULL)
  
  return(IEs)
}

## (4.3) plot the IE trajectories
case1_df_IEs_for_plotting <- apply(subsamples, 1, case1_get_IEs_for_plotting) 
case4_df_IEs_for_plotting <- apply(subsamples, 1, case4_get_IEs_for_plotting) 

case1_df_IEs_for_plotting <- do.call(rbind, case1_df_IEs_for_plotting)  
case4_df_IEs_for_plotting <- do.call(rbind, case4_df_IEs_for_plotting)  

case1_subsamples_IE_infection <- ggplot() +
  geom_line(data = case1_df_IEs_for_plotting, aes(x = t, y= IE_infection, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.01,0.05), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = "none")

case1_subsamples_IE_death <- ggplot() +
  geom_line(data = case1_df_IEs_for_plotting, aes(x = t, y= IE_death, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.0001,0.0005), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = "none")

case4_subsamples_IE_infection <- ggplot() +
  geom_line(data = case4_df_IEs_for_plotting, aes(x = t, y= IE_infection, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.01,0.05), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = "none")

case4_subsamples_IE_death <- ggplot() +
  geom_line(data = case4_df_IEs_for_plotting, aes(x = t, y= IE_death, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.0001,0.0005), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = "none")

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 4", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot3 <- "
#fffffggggg
ahhhhhiiiii
ahhhhhiiiii
ahhhhhiiiii
bjjjjjkkkkk
bjjjjjkkkkk
bjjjjjkkkkk
"

plotlist3 <-
  list(
    a = row1,
    b = row2,
    f = col1,
    g = col2,
    h = case1_subsamples_IE_infection,
    i = case1_subsamples_IE_death,
    j = case4_subsamples_IE_infection,
    k = case4_subsamples_IE_death)

plot_subsamples_IE_labeled <- wrap_plots(plotlist3, guides = 'collect', design = layoutplot3) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig2.png", case1_subsamples_IE_death, width = 10, height=5, dpi=300, units="in")
