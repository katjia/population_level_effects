# <I> definitions --------------------------------------------------------------
## alpha_observed = "observed" proportion vaccinated
## alpha_hypothetical = "hypothetical" alpha

# <II> functions ---------------------------------------------------------------
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
} # checked

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
} # checked

## a function to left join PDE with POE
get_PDE_POE <- function(dataframe, alpha_hypothetical, alpha_observed){
  
  dataframe_PDE <- get_PDE(dataframe, alpha_hypothetical, alpha_observed)

  dataframe_POE <- get_POE(dataframe)

  dataframe_PDE %>%
    select(t, PDE_infection, PDE_death) %>%
    left_join(dataframe_POE, by = "t")
} # checked

# <III> verify Claim 1a (alpha_hypothetical = 0, alpha_observed = 0.7) ---------
## (1) get PDE and POE ---------------------------------------------------------
PDE_POE_alpha0_0_alpha1_07 <-
  lapply(df_all_counterfactuals_scenarios_1_to_5, function(df) {
    df %>%
      filter(alpha %in% c(0, 0.7)) %>%
      get_PDE_POE(alpha_hypothetical = 0, alpha_observed = 0.7)
  })

## (2) plot PDE and POE --------------------------------------------------------
plot_claim1a_PDE_POE_infection <- lapply(PDE_POE_alpha0_0_alpha1_07, function(df) {
    ggplot(df) +
      geom_line(aes(x = t, y= POE_infection, lty = "Overall")) +
      geom_line(aes(x = t, y= PDE_infection, lty = "Direct")) +
      labs(y=bquote(atop('POE'^'infection'*'('*t*','*0*','*0.7*')'~or,
                         ~'PDE'^'infection'*'('*t*','*0*','*0.7*')')),
           x="Day",
           lty=" ") +
      theme_bw() +
      theme(legend.position = "none")
  })

plot_claim1a_PDE_POE_death <- lapply(PDE_POE_alpha0_0_alpha1_07, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= POE_death, lty = "Overall")) +
    geom_line(aes(x = t, y= PDE_death, lty = "Direct")) +
    labs(y=bquote(atop('POE'^'death'*'('*t*','*0*','*0.7*')'~or,
                       ~'PDE'^'death'*'('*t*','*0*','*0.7*')')),
         x="Day",
         lty=" ") +
    theme_bw() +
    theme(legend.direction = "horizontal",
          legend.text=element_text(size=12))
})

## (3) collate plots -----------------------------------------------------------
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 5", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_fig2 <- "
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

plotlist_fig2 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_claim1a_PDE_POE_infection[["scenario_1"]],
    i = plot_claim1a_PDE_POE_death[["scenario_1"]] + labs(caption = NULL),
    j = plot_claim1a_PDE_POE_infection[["scenario_2"]],
    k = plot_claim1a_PDE_POE_death[["scenario_2"]] + labs(caption = NULL),
    l = plot_claim1a_PDE_POE_infection[["scenario_3"]],
    m = plot_claim1a_PDE_POE_death[["scenario_3"]] + labs(caption = NULL),
    n = plot_claim1a_PDE_POE_infection[["scenario_4"]],
    o = plot_claim1a_PDE_POE_death[["scenario_4"]] + labs(caption = NULL),
    p = plot_claim1a_PDE_POE_infection[["scenario_5"]],
    q = plot_claim1a_PDE_POE_death[["scenario_5"]] + labs(caption = NULL),
    z = guide_area())

plot_claim1a_PDE_POE_labeled <- wrap_plots(plotlist_fig2, guides = 'collect', design = layoutplot_fig2) 

ggsave("3_figures/Fig2.png", plot_claim1a_PDE_POE_labeled, width = 9, height=10, dpi=300, units="in")

# <IV> Verify Claim 1b (alpha = 0.9, alpha_observed = 0.7) ---------------------
## (1) get PDE and POE ---------------------------------------------------------
PDE_POE_alpha1_07_alpha2_09 <-
  lapply(df_all_counterfactuals_scenarios_1_to_5, function(df) {
    df %>%
      filter(alpha %in% c(0.9, 0.7)) %>%
      get_PDE_POE(alpha_hypothetical = 0.9, alpha_observed = 0.7)
  })

## (2) plot PDE and POE --------------------------------------------------------
plot_claim1b_PDE_POE_infection <- lapply(PDE_POE_alpha1_07_alpha2_09, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= POE_infection, lty = "Overall")) +
    geom_line(aes(x = t, y= PDE_infection, lty = "Direct")) +
    labs(y=bquote(atop('POE'^'infection'*'('*t*','*0.7*','*0.9*')'~or,
                       ~'PDE'^'infection'*'('*t*','*0.9*','*0.7*')')),
         x="Day",
         lty=" ") +
    theme_bw() +
    theme(legend.position = "none")
})

plot_claim1b_PDE_POE_death <- lapply(PDE_POE_alpha1_07_alpha2_09, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= POE_death, lty = "Overall")) +
    geom_line(aes(x = t, y= PDE_death, lty = "Direct")) +
    labs(y=bquote(atop('POE'^'death'*'('*t*','*0.7*','*0.9*')'~or,
                       ~'PDE'^'death'*'('*t*','*0.9*','*0.7*')')),
         x="Day",
         lty=" ") +
    theme_bw() +
    theme(legend.direction = "horizontal",
          legend.text=element_text(size=12))
})

## (3) collate plots  -----------------------------------------------------------
plotlist_fig3 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_claim1b_PDE_POE_infection[["scenario_1"]],
    i = plot_claim1b_PDE_POE_death[["scenario_1"]] + labs(caption = NULL),
    j = plot_claim1b_PDE_POE_infection[["scenario_2"]],
    k = plot_claim1b_PDE_POE_death[["scenario_2"]] + labs(caption = NULL),
    l = plot_claim1b_PDE_POE_infection[["scenario_3"]],
    m = plot_claim1b_PDE_POE_death[["scenario_3"]] + labs(caption = NULL),
    n = plot_claim1b_PDE_POE_infection[["scenario_4"]],
    o = plot_claim1b_PDE_POE_death[["scenario_4"]] + labs(caption = NULL),
    p = plot_claim1b_PDE_POE_infection[["scenario_5"]],
    q = plot_claim1b_PDE_POE_death[["scenario_5"]] + labs(caption = NULL),
    z = guide_area())

plot_claim1b_PDE_POE_labeled <- wrap_plots(plotlist_fig3, guides = 'collect', design = layoutplot_fig2) 

ggsave("3_figures/Fig3.png", plot_claim1b_PDE_POE_labeled, width = 9, height=10, dpi=300, units="in")
