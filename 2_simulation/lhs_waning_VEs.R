# latin hypercube sampling to test whether Claim 1 holds under Case 1 and Case 4
library(lhs)
library(ggforce)
set.seed(123)

# (1) LHS ----------------------------------------------------------------------
ranges <- list(
  "R0" = c(1, 20),
  "alpha" = c(0.01, 1),
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
parms_values <- as.data.frame(randomLHS(ndraws, 10)) 
names(parms_values) <- names(ranges)
for(parm in names(ranges)){
  parms_values[,parm] <- parms_values[,parm]*(ranges[[parm]][2]-ranges[[parm]][1]) + ranges[[parm]][1]
}

# (2) Test Claim 1 under Case 1 ------------------------------------------------

# call Odin models
case1_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case1_time_invariant_par.R")
case4_generator <- odin::odin("~/Documents/GitHub/population_level_effects/1_model/case4_waning_VEs.R")

# a function to get pInf_s and pDeath_s for any proportion vaccinated
get_pInf_s_pDeath_s_any_alpha <- function(dataframe){
  dataframe %>%
    select(alpha, t, Cum_inf_s, D_s) %>%
    mutate(
      pInf_s = case_when(alpha > 0 ~ Cum_inf_s / (2e4 * (1 - alpha)),
                         alpha == 0.0 ~ Cum_inf_s / 2e4),
      pDeath_s = case_when(alpha > 0 ~ D_s / (2e4 * (1 - alpha)),
                           alpha == 0.0 ~ D_s / 2e4))}

# a function to get IE infection for any proportion vaccinated
get_IE_infection_any_alpha <- function(dataframe) {
  dataframe %>%
    select(alpha, t, pInf_s) %>%
    group_by(t) %>%
    summarise(IE_infection = last(pInf_s) - first(pInf_s))
}

# a function to get IE death for any proportion vaccinated
get_IE_death_any_alpha <- function(dataframe) {
  dataframe %>%
    select(alpha, t, pDeath_s) %>%
    group_by(t) %>%
    summarise(IE_death = last(pDeath_s) - first(pDeath_s))
}

# a function to run model under case 1
simulate_case1 <- function(par){
  
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
                                              alpha = alpha))}, .id = "par")
  return(out)}

# a function to test Claim 1 under case 1
test_claim1_case1 <- function(par){
  
  # simulate
  out <- simulate_case1(par)
  
  # get pInf_s and pDeath_s
  pInf_s_pDeath_s <- get_pInf_s_pDeath_s_any_alpha(out) 

  # get IE infection 
  IE_infection <- get_IE_infection_any_alpha(pInf_s_pDeath_s)
  
  days_neg_IE_infection <- sum(IE_infection$IE_infection < -10e-7)
  
  # get IE death
  IE_death <- get_IE_death_any_alpha(pInf_s_pDeath_s)
  
  days_neg_IE_death <- sum(IE_death$IE_death < -10e-7)
  
  return(c(days_neg_IE_infection=days_neg_IE_infection,
           days_neg_IE_death=days_neg_IE_death))
}

## (2.1) Total days with neg IE
case1_days_neg_IE <- apply(parms_values, 1, test_claim1_case1)
colnames(case1_days_neg_IE) <- 1:ndraws

## (2.2) Identify par combo for which IE is negative
(case1_neg_IE_infection_par <- which(case1_days_neg_IE[1,] > 0, TRUE))
(case1_neg_IE_death_par <- which(case1_days_neg_IE[2,] > 0, TRUE))

# (3) Test Claim 1 under Case 4 ------------------------------------------------
# a function to run model under case 4
simulate_case4 <- function(par){

  # Take par values from LHS
  prop_vax <- as.numeric(par["alpha"])
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
  
  # Simulate counterfactuals
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
  
  return(out)
}

# a function to test Claim 1 under case 4
test_claim1_case4 <- function(par){
  
  # simulate
  out <- simulate_case4(par)
  
  # a function to get pInf_s and pDeath_s
  pInf_s_pDeath_s <- get_pInf_s_pDeath_s_any_alpha(out) 
  
  # get IE infection 
  IE_infection <- get_IE_infection_any_alpha(pInf_s_pDeath_s)
  
  days_neg_IE_infection <- sum(IE_infection$IE_infection < -10e-7)
  
  # get IE death
  IE_death <- get_IE_death_any_alpha(pInf_s_pDeath_s)
  
  days_neg_IE_death <- sum(IE_death$IE_death < -10e-7)
  
  return(c(days_neg_IE_infection=days_neg_IE_infection,
           days_neg_IE_death=days_neg_IE_death))
}

## (3.1) Total days with neg IE
case4_days_neg_IE <- apply(parms_values, 1, test_claim1_case4)
colnames(case4_days_neg_IE) <- 1:ndraws

## (3.2) Identify par combo for which IE is negative
(case4_neg_IE_infection_par <- which(case4_days_neg_IE[1,] > 0, TRUE))
(case4_neg_IE_death_par <- which(case4_days_neg_IE[2,] > 0, TRUE))

# (4) Randomly draw n samples from LHS and plot the IEs ------------------------
## (4.1) randomly draw samples
set.seed(123)
subsamples_id <- sample(1:ndraws, size=50, replace = FALSE)
subsamples <- parms_values[subsamples_id, ]
subsamples <- cbind(subsamples, par_id=subsamples_id)

## (4.2) get IEs for plotting
# a function to get IE for plotting under Case 1
get_IEs_for_plotting_case1 <- function(subpar){
  
  out <- simulate_case1(subpar)
  
  # get pInf_s and pDeath_s
  pInf_s_pDeath_s <- get_pInf_s_pDeath_s_any_alpha(out) 
  
  # get IE infection 
  IE_infection <- get_IE_infection_any_alpha(pInf_s_pDeath_s)
  
  # get IE death
  IE_death <- get_IE_death_any_alpha(pInf_s_pDeath_s)
  
  IEs <- IE_death %>%
    left_join(IE_infection, by="t") %>%
    cbind(alpha1=prop_vax,
          par_id=subpar["par_id"],
          row.names = NULL)
  
  return(IEs)
  
}

# a function to get IE for plotting under Case 4
get_IEs_for_plotting_case4 <- function(subpar){
  
  out <- simulate_case4(subpar)
  
  # get pInf_s and pDeath_s
  pInf_s_pDeath_s <- get_pInf_s_pDeath_s_any_alpha(out) 
  
  # get IE infection 
  IE_infection <- get_IE_infection_any_alpha(pInf_s_pDeath_s)
  
  # get IE death
  IE_death <- get_IE_death_any_alpha(pInf_s_pDeath_s)
  
  IEs <- IE_death %>%
    left_join(IE_infection, by="t") %>%
    cbind(alpha1=prop_vax,
          par_id=subpar["par_id"],
          row.names = NULL)
  
  return(IEs)
  
}

## (4.3) plot the IE trajectories
df_IEs_for_plotting_case1 <- apply(subsamples, 1, get_IEs_for_plotting_case1) 
df_IEs_for_plotting_case4 <- apply(subsamples, 1, get_IEs_for_plotting_case4) 

df_IEs_for_plotting_case1 <- do.call(rbind, df_IEs_for_plotting_case1)  
df_IEs_for_plotting_case4 <- do.call(rbind, df_IEs_for_plotting_case4)  

subsamples_IE_infection_case1 <- ggplot() +
  geom_line(data = df_IEs_for_plotting_case1, aes(x = t, y= IE_infection, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'infection'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(zoom.size = 1, ylim = c(-0.01,0.05)) +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

subsamples_IE_death_case1 <- ggplot() +
  geom_line(data = df_IEs_for_plotting_case1, aes(x = t, y= IE_death, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'death'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(zoom.size = 1, ylim = c(-0.0001,0.0005)) +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

subsamples_IE_infection_case4 <- ggplot() +
  geom_line(data = df_IEs_for_plotting_case4, aes(x = t, y= IE_infection, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'infection'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(zoom.size = 1, ylim = c(-0.01,0.05)) +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

subsamples_IE_death_case4 <- ggplot() +
  geom_line(data = df_IEs_for_plotting_case4, aes(x = t, y= IE_death, group = par_id), col="gray40", alpha=0.5, linewidth=0.25) +
  labs(tag=" ",
       y=bquote('IIE'^'death'*'('*t*','*0*','*alpha[1]*')'),
       x="Day") +
  facet_zoom(zoom.size = 1, ylim = c(-0.0001,0.0005)) +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 4", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot3 <- "
#fffffffggggggg
ahhhhhhhiiiiiii
ahhhhhhhiiiiiii
ahhhhhhhiiiiiii
bjjjjjjjkkkkkkk
bjjjjjjjkkkkkkk
bjjjjjjjkkkkkkk
"

plotlist3 <-
  list(
    a = row1,
    b = row2,
    f = col1,
    g = col2,
    h = subsamples_IE_infection_case1,
    i = subsamples_IE_death_case1,
    j = subsamples_IE_infection_case4,
    k = subsamples_IE_death_case4)

plot_subsamples_IE_labeled <- wrap_plots(plotlist3, guides = 'collect', design = layoutplot3) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig2_2.png", plot_subsamples_IE_labeled, width = 11, height=5, dpi=300, units="in")
