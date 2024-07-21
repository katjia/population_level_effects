library(ggforce)
# <I> DEFINITIONS ---------------------------------------------------------------
## alpha_observed = "observed" proportion vaccinated
## alpha_hypothetical = "hypothetical" alpha

## a function to get PDE
get_PDE <- function(dataframe, alpha_hypothetical, alpha_observed) {
  abs_alpha_diff <- abs(alpha_hypothetical - alpha_observed)
  dataframe %>%
    filter(alpha == alpha_observed) %>%
    mutate(
      pInf_s = Cum_inf_s / (2e4 * (1 - alpha)),
      pInf_v = Cum_inf_v / (2e4 * alpha),
      pDeath_s = D_s / (2e4 * (1 - alpha)),
      pDeath_v = D_v / (2e4 * alpha),
      PDE_death = (2e4 * abs_alpha_diff) * (pDeath_s - pDeath_v),
      PDE_infection = (2e4 * abs_alpha_diff) * (pInf_s - pInf_v)
    )
}

# <II> FUNCTIONS ---------------------------------------------------------------
## a function to get POE
get_POE <- function(dataframe) {
  ### 1. get overall prob
  df_pInf_pDeath <- dataframe %>%
    mutate(pDeath = D / 2e4,
           pInf = Cum_inf / 2e4)
  ### 2. get POE
  df_pInf_pDeath %>%
    select(alpha, t, pDeath, pInf) %>%
    group_by(t) %>%
    arrange(t,-alpha) %>%
    summarise(POE_infection = 2e4 * (last(pInf) - first(pInf)),
              POE_death = 2e4 * (last(pDeath) - first(pDeath)))
}

# <III> Verify Claim 1a (alpha = 0, alpha_observed = 0.7) --------------------------
## (1) Case 1: Time-invariant parameters ---------------------------------------
### (1.1) get PDE
df_alpha0_0_alpha1_07_PDE_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0, alpha_observed = 0.7)

## (1.2) get POE
df_alpha0_0_alpha1_07_POE_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_POE()

## (1.3) plot
plot_claim1a_PDE_POE_infection_case1 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case1, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case1, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1a_PDE_POE_death_case1 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case1, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case1, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (2) Case 2: Increasing effective contacts (beta) -----------------------------
### (2.1) get PDE
df_alpha0_0_alpha1_07_PDE_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0, alpha_observed = 0.7)

### (2.2) get POE 
df_alpha0_0_alpha1_07_POE_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_POE()

### (2.3) plot
plot_claim1a_PDE_POE_infection_case2 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case2, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case2, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1a_PDE_POE_death_case2 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case2, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case2, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (3) Case 3: Increasing mu ----------------------------------------------------
### (3.1) get PDE
df_alpha0_0_alpha1_07_PDE_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0, alpha_observed = 0.7)

### (3.2) get POE 
df_alpha0_0_alpha1_07_POE_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_POE()

### (3.3) plot 
plot_claim1a_PDE_POE_infection_case3 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case3, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case3, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1a_PDE_POE_death_case3 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case3, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case3, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (4) Case 4: Waning VEs -------------------------------------------------------
### (4.1) get PDE
df_alpha0_0_alpha1_07_PDE_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0, alpha_observed = 0.7)

### (4.2) get POE 
df_alpha0_0_alpha1_07_POE_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_POE()

### (4.3) plot
plot_claim1a_PDE_POE_infection_case4 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case4, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case4, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1a_PDE_POE_death_case4 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case4, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case4, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (5) Case 5: Waning + increasing beta --------------------------------------------
### (5.1) get PDE
df_alpha0_0_alpha1_07_PDE_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0, alpha_observed = 0.7)

### (5.2) get POE 
df_alpha0_0_alpha1_07_POE_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_POE()

### (5.3) plot 
plot_claim1a_PDE_POE_infection_case5 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case5, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case5, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,~
                     'PDE'^'infection'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1a_PDE_POE_death_case5 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_POE_case5, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha0_0_alpha1_07_PDE_case5, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                     ~'PDE'^'death'*'('*t*','*0*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# <IV> Plot Figure 1 -----------------------------------------------------------
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 5", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_fig1 <- "
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

plotlist_fig1 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_claim1a_PDE_POE_infection_case1,
    i = plot_claim1a_PDE_POE_death_case1 + labs(caption = NULL),
    j = plot_claim1a_PDE_POE_infection_case2,
    k = plot_claim1a_PDE_POE_death_case2 + labs(caption = NULL),
    l = plot_claim1a_PDE_POE_infection_case3,
    m = plot_claim1a_PDE_POE_death_case3 + labs(caption = NULL),
    n = plot_claim1a_PDE_POE_infection_case4,
    o = plot_claim1a_PDE_POE_death_case4 + labs(caption = NULL),
    p = plot_claim1a_PDE_POE_infection_case5,
    q = plot_claim1a_PDE_POE_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_claim1a_PDE_POE_labeled <- wrap_plots(plotlist_fig1, guides = 'collect', design = layoutplot_fig1) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/Fig1.png", plot_claim1a_PDE_POE_labeled, width = 9, height=10, dpi=300, units="in")

# <V> Verify Claim 1b (alpha = 0.9, alpha_observed = 0.7) --------------------------
## (1) Case 1: Time-invariant parameters ---------------------------------------
### (1.1) get PDE
df_alpha1_07_alpha2_09_PDE_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0.9, alpha_observed = 0.7)

## (1.2) get POE
df_alpha1_07_alpha2_09_POE_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_POE()

## (1.3) plot
plot_claim1b_PDE_POE_infection_case1 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case1, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case1, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1b_PDE_POE_death_case1 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case1, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case1, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (2) Case 2: Increasing effective contacts (beta) ----------------------------
### (2.1) get PDE
df_alpha1_07_alpha2_09_PDE_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0.9, alpha_observed = 0.7)

### (2.2) get POE 
df_alpha1_07_alpha2_09_POE_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_POE()

### (2.3) plot
plot_claim1b_PDE_POE_infection_case2 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case2, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case2, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1b_PDE_POE_death_case2 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case2, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case2, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (3) Case 3: Increasing mu ----------------------------------------------------
### (3.1) get PDE
df_alpha1_07_alpha2_09_PDE_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_PDE(alpha_hypothetical = 0.9, alpha_observed = 0.7)

### (3.2) get POE 
df_alpha1_07_alpha2_09_POE_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_POE()

### (3.3) plot 
plot_claim1b_PDE_POE_infection_case3 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case3, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case3, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1b_PDE_POE_death_case3 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case3, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case3, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (4) Case 4: Waning VEs -------------------------------------------------------
### (4.1) get PDE
df_alpha1_07_alpha2_09_PDE_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_PDE(alpha = 0.9, alpha_observed = 0.7)

### (4.2) get POE 
df_alpha1_07_alpha2_09_POE_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_POE()

### (4.3) plot
plot_claim1b_PDE_POE_infection_case4 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case4, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case4, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1b_PDE_POE_death_case4 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case4, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case4, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

## (5) Case 5: Waning + increasing beta --------------------------------------------
### (5.1) get PDE
df_alpha1_07_alpha2_09_PDE_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_PDE(alpha = 0.9, alpha_observed = 0.7)

### (5.2) get POE 
df_alpha1_07_alpha2_09_POE_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0.9, 0.7)) %>%
  get_POE()

### (5.3) plot 
plot_claim1b_PDE_POE_infection_case5 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case5, aes(x = t, y= POE_infection, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case5, aes(x = t, y= PDE_infection, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'infection'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'infection'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.position = "none")

plot_claim1b_PDE_POE_death_case5 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_POE_case5, aes(x = t, y= POE_death, lty = "Overall")) +
  geom_line(data = df_alpha1_07_alpha2_09_PDE_case5, aes(x = t, y= PDE_death, lty = "Direct")) +
  labs(y=bquote(atop('POE'^'death'*'('*t*','*0.7*','*0.9*')'~or,
                     ~'PDE'^'death'*'('*t*','*0.9*','*0.7*')')),
       x="Day",
       lty=" ") +
  theme_bw() +
  theme(legend.direction = 'horizontal')

# <VI> Plot Figure 2 -----------------------------------------------------------
plotlist_fig2 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_claim1b_PDE_POE_infection_case1,
    i = plot_claim1b_PDE_POE_death_case1 + labs(caption = NULL),
    j = plot_claim1b_PDE_POE_infection_case2,
    k = plot_claim1b_PDE_POE_death_case2 + labs(caption = NULL),
    l = plot_claim1b_PDE_POE_infection_case3,
    m = plot_claim1b_PDE_POE_death_case3 + labs(caption = NULL),
    n = plot_claim1b_PDE_POE_infection_case4,
    o = plot_claim1b_PDE_POE_death_case4 + labs(caption = NULL),
    p = plot_claim1b_PDE_POE_infection_case5,
    q = plot_claim1b_PDE_POE_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_claim1b_PDE_POE_labeled <- wrap_plots(plotlist_fig2, guides = 'collect', design = layoutplot_fig1) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/Fig2.png", plot_claim1b_PDE_POE_labeled, width = 9, height=10, dpi=300, units="in")
