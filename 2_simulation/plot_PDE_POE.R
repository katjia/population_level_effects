library(ggforce)
# (1) Case 1: Time-invariant parameters ----------------------------------------
## (1.1) get PDE
### a function to get PDE
get_PDE <- function(dataframe){
  dataframe %>%
    filter(alpha %in% 0.7) %>%
    filter(VE_infection %in% c(0.5)) %>%
    mutate(
      pInf_s = Cum_inf_s / (2e4 * 0.3),
      pInf_v = Cum_inf_v / (2e4 * 0.7),
      pDeath_s = D_s / (2e4 * 0.3),
      pDeath_v = D_v / (2e4 * 0.7),
      PDE_death = (2e4 * 0.7) * (pDeath_s - pDeath_v),
      PDE_infection = (2e4 * 0.7) * (pInf_s - pInf_v)
    )
}

df_PDE_case1 <- get_PDE(df_all_counterfactuals_case1_time_invariant_par)

## (1.2) get POE 
### a function to get POE
get_POE <- function(dataframe){
  ### 1. get overall prob
  df_pInf_pDeath <- dataframe %>%
    filter(VE_infection %in% c(0.5)) %>%
    mutate(pDeath=D/2e4,
           pInf=Cum_inf/2e4)
  ### 2. get POE 
  df_pInf_pDeath %>%
    select(VE_infection, alpha, t, pDeath, pInf) %>%
    group_by(VE_infection, t) %>%
    pivot_wider(names_from = alpha, values_from = c(pDeath, pInf)) %>%
    mutate(POE_death = 2e4 * (pDeath_0 - pDeath_0.7),
           POE_infection = 2e4 * (pInf_0 - pInf_0.7)) 
}

df_POE_case1 <- get_POE(df_all_counterfactuals_case1_time_invariant_par)

## (1.3) plot
plot_PDE_POE_infection_case1 <- ggplot() +
  geom_line(data = df_POE_case1, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_PDE_case1, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(tag="A(i)",
       y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_PDE_POE_death_case1 <- ggplot() +
  geom_line(data = df_POE_case1, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_PDE_case1, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(tag="A(ii)",
       y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case1_time_invariant_par.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_PDE_POE.R")) +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# (2) Case 2: Increasing effective contacts (beta) -----------------------------
## (2.1) get PDE
df_PDE_case2 <- get_PDE(df_all_counterfactuals_case2_inc_beta)

## (2.2) get POE 
df_POE_case2 <- get_POE(df_all_counterfactuals_case2_inc_beta)

## (2.3) plot
plot_PDE_POE_infection_case2 <- ggplot() +
  geom_line(data = df_POE_case2, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_PDE_case2, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(tag="B(i)",
       y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_PDE_POE_death_case2 <- ggplot() +
  geom_line(data = df_POE_case2, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_PDE_case2, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(tag="B(ii)",
       y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ",
       caption = paste0("~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/1_model/case2_inc_beta.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_PDE_POE.R")) +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# (3) Case 3: Increasing mu ----------------------------------------------------
## (3.1) get PDE
df_PDE_case3 <- get_PDE(df_all_counterfactuals_case3_inc_mu)

## (3.2) get POE 
### (3.2.1) compute overall prob
df_POE_case3 <- get_POE(df_all_counterfactuals_case3_inc_mu)

## (3.3) plot 
plot_PDE_POE_infection_case3 <- ggplot() +
  geom_line(data = df_POE_case3, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_PDE_case3, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(tag="C(i)",
       y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_PDE_POE_death_case3 <- ggplot() +
  geom_line(data = df_POE_case3, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_PDE_case3, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(tag="C(ii)",
       y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case3_inc_mu.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_PDE_POE.R")) +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# (4) Case 4: Waning VEs --------------------------------------------
## (4.1) get PDE
df_PDE_case4 <- get_PDE(df_all_counterfactuals_case4_waning)

## (4.2) get POE 
df_POE_case4 <- get_POE(df_all_counterfactuals_case4_waning)

## (4.3) plot
plot_PDE_POE_infection_case4 <- ggplot() +
  geom_line(data = df_POE_case4, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_PDE_case4, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(tag="D(i)",
       y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_PDE_POE_death_case4 <- ggplot() +
  geom_line(data = df_POE_case4, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_PDE_case4, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(tag="D(ii)", 
       y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case4_waning_VEs.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_PDE_POE.R")) +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# Case 5: Waning + increasing beta --------------------------------------------
## (5.1) get PDE
df_PDE_case5 <- get_PDE(df_all_counterfactuals_case5_inc_beta_waning)

## (5.2) get POE 
df_POE_case5 <- get_POE(df_all_counterfactuals_case5_inc_beta_waning)

## (5.3) plot 
plot_PDE_POE_infection_case5 <- ggplot() +
  geom_line(data = df_POE_case5, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_PDE_case5, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(tag="E(i)",
       y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,~
                     'PDE'^'infection'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_PDE_POE_death_case5 <- ggplot() +
  geom_line(data = df_POE_case5, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_PDE_case5, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(tag="E(ii)",
       y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.7*')')),
       x="Day",
       lty=" ",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case5_inc_beta_waning_VEs.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_PDE_POE.R")) +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# collate the plots
layoutplot2 <- "
#ffffffgggggg
ahhhhhhiiiiii
ahhhhhhiiiiii
ahhhhhhiiiiii
bjjjjjjkkkkkk
bjjjjjjkkkkkk
bjjjjjjkkkkkk
cllllllmmmmmm
cllllllmmmmmm
cllllllmmmmmm
dnnnnnnoooooo
dnnnnnnoooooo
dnnnnnnoooooo
eppppppqqqqqq
eppppppqqqqqq
eppppppqqqqqq
####zzzzz####
"

plotlist2 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_PDE_POE_infection_case1,
    i = plot_PDE_POE_death_case1 + labs(caption = NULL),
    j = plot_PDE_POE_infection_case2,
    k = plot_PDE_POE_death_case2 + labs(caption = NULL),
    l = plot_PDE_POE_infection_case3,
    m = plot_PDE_POE_death_case3 + labs(caption = NULL),
    n = plot_PDE_POE_infection_case4,
    o = plot_PDE_POE_death_case4 + labs(caption = NULL),
    p = plot_PDE_POE_infection_case5,
    q = plot_PDE_POE_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_PDE_POE_labeled <- wrap_plots(plotlist2, guides = 'collect', design = layoutplot2) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/Fig2_PDE_POE_2.png", plot_PDE_POE_labeled, width = 9, height=10, dpi=300, units="in")
