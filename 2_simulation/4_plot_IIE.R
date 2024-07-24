# <I> functions ----------------------------------------------------------------
### create a function to get IIE^Unvax
get_IIE_unvax <- function(dataframe){
  dataframe %>%
    select(-par, alpha, t, Cum_inf_s, D_s) %>%
    mutate(pInf_s = Cum_inf_s / (2e4 * (1-alpha)), 
           pDeath_s = D_s / (2e4 * (1-alpha))) %>%
    select(alpha, t, pInf_s, pDeath_s) %>%
    group_by(t) %>%
    arrange(t,-alpha) %>%
    summarise(IIE_infection = last(pInf_s) - first(pInf_s),
              IIE_death = last(pDeath_s) - first(pDeath_s))
}

### create a function to get IIE^Vax
get_IIE_vax <- function(dataframe){
  dataframe %>%
    select(-par, alpha, t, Cum_inf_v, D_v) %>%
    mutate(pInf_v = Cum_inf_v / (2e4 * alpha), 
           pDeath_v = D_v / (2e4 * alpha)) %>%
    select(alpha, t, pInf_v, pDeath_v) %>%
    group_by(t) %>%
    arrange(t,-alpha) %>%
    summarise(IIE_infection = last(pInf_v) - first(pInf_v),
              IIE_death = last(pDeath_v) - first(pDeath_v))
}

# <II> plot IIE^Unvax(t, alpha0=0, alpha1=0.7) for Claim 1a --------------------
## (1) get IIE unvax -----------------------------------------------------------
IIE_unvax_alpha0_0_alpha1_07 <-
  lapply(df_all_counterfactuals_scenarios_1_to_5, function(df) {
    df %>%
      filter(alpha %in% c(0, 0.7)) %>%
      get_IIE_unvax()
  })
## (2) plot IIE unvax ----------------------------------------------------------
plot_alpha0_0_alpha1_07_IIE_unvax_infection <- lapply(IIE_unvax_alpha0_0_alpha1_07, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= IIE_infection)) +
    labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0*','*0.7*')'),
         x="Day") +
    theme_bw() +
    theme(legend.position = "none",
          zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
          zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 
})

plot_alpha0_0_alpha1_07_IIE_unvax_death <- lapply(IIE_unvax_alpha0_0_alpha1_07, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= IIE_death)) +
    labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0*','*0.7*')'),
         x="Day") +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal',
          zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
          zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 
})

## (3) collate the plots -------------------------------------------------------
layoutplot_eFig3 <- "
#fffffffggggggg
ahhhhhhhiiiiiii
ahhhhhhhiiiiiii
bjjjjjjjkkkkkkk
bjjjjjjjkkkkkkk
clllllllmmmmmmm
clllllllmmmmmmm
dnnnnnnnooooooo
dnnnnnnnooooooo
epppppppqqqqqqq
epppppppqqqqqqq
#####zzzzz#####
"

plotlist_eFig3 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_alpha0_0_alpha1_07_IIE_unvax_infection$scenario_1,
    i = plot_alpha0_0_alpha1_07_IIE_unvax_death$scenario_1 + labs(caption = NULL),
    j = plot_alpha0_0_alpha1_07_IIE_unvax_infection$scenario_2,
    k = plot_alpha0_0_alpha1_07_IIE_unvax_death$scenario_2 + labs(caption = NULL),
    l = plot_alpha0_0_alpha1_07_IIE_unvax_infection$scenario_3,
    m = plot_alpha0_0_alpha1_07_IIE_unvax_death$scenario_3 + labs(caption = NULL),
    n = plot_alpha0_0_alpha1_07_IIE_unvax_infection$scenario_4,
    o = plot_alpha0_0_alpha1_07_IIE_unvax_death$scenario_4 + labs(caption = NULL),
    p = plot_alpha0_0_alpha1_07_IIE_unvax_infection$scenario_5,
    q = plot_alpha0_0_alpha1_07_IIE_unvax_death$scenario_5 + labs(caption = NULL),
    z = guide_area())

plot_alpha0_0_alpha1_07_IIE_unvax_labeled <- wrap_plots(plotlist_eFig3, guides = 'collect', design = layoutplot_eFig3) 

ggsave("3_figures/eFig3.png", plot_alpha0_0_alpha1_07_IIE_unvax_labeled, width = 12, height=12, dpi=300, units="in")

# <III> plot IIE^Unvax(t, alpha1=0.7, alpha2=0.9) for Claim 1b -----------------
## (1) get IIE unvax -----------------------------------------------------------
IIE_unvax_alpha1_07_alpha2_09 <-
  lapply(df_all_counterfactuals_scenarios_1_to_5, function(df) {
    df %>%
      filter(alpha %in% c(0.7, 0.9)) %>%
      get_IIE_unvax()
  })

## (2) plot IIE unvax ----------------------------------------------------------
plot_alpha1_07_alpha2_09_IIE_unvax_infection <- lapply(IIE_unvax_alpha1_07_alpha2_09, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= IIE_infection)) +
    labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0.7*','*0.9*')'),
         x="Day") +
    theme_bw() +
    theme(legend.position = "none",
          zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
          zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 
})

plot_alpha1_07_alpha2_09_IIE_unvax_death <- lapply(IIE_unvax_alpha1_07_alpha2_09, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= IIE_death)) +
    labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0.7*','*0.9*')'),
         x="Day") +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal',
          zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
          zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 
})

## (3) collate the plots -------------------------------------------------------
plotlist_eFig4 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_alpha1_07_alpha2_09_IIE_unvax_infection$scenario_1,
    i = plot_alpha1_07_alpha2_09_IIE_unvax_death$scenario_1 + labs(caption = NULL),
    j = plot_alpha1_07_alpha2_09_IIE_unvax_infection$scenario_2,
    k = plot_alpha1_07_alpha2_09_IIE_unvax_death$scenario_2 + labs(caption = NULL),
    l = plot_alpha1_07_alpha2_09_IIE_unvax_infection$scenario_3,
    m = plot_alpha1_07_alpha2_09_IIE_unvax_death$scenario_3 + labs(caption = NULL),
    n = plot_alpha1_07_alpha2_09_IIE_unvax_infection$scenario_4,
    o = plot_alpha1_07_alpha2_09_IIE_unvax_death$scenario_4 + labs(caption = NULL),
    p = plot_alpha1_07_alpha2_09_IIE_unvax_infection$scenario_5,
    q = plot_alpha1_07_alpha2_09_IIE_unvax_death$scenario_5 + labs(caption = NULL),
    z = guide_area())

plot_alpha1_07_alpha2_09_IIE_unvax_labeled <- wrap_plots(plotlist_eFig4, guides = 'collect', design = layoutplot_eFig3) 

ggsave("3_figures/eFig4.png", plot_alpha1_07_alpha2_09_IIE_unvax_labeled, width = 12, height=12, dpi=300, units="in")

# <IV> plot IIE^Vax(t, alpha1=0.7, alpha2=0.9) for Claim 1b --------------------
## (1) get IIE vax -------------------------------------------------------------
IIE_vax_alpha1_07_alpha2_09 <-
  lapply(df_all_counterfactuals_scenarios_1_to_5, function(df) {
    df %>%
      filter(alpha %in% c(0.7, 0.9)) %>%
      get_IIE_vax()
  })

## (2) plot IIE vax ----------------------------------------------------------
plot_alpha1_07_alpha2_09_IIE_vax_infection <- lapply(IIE_vax_alpha1_07_alpha2_09, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= IIE_infection)) +
    labs(y=bquote('IIE'^'Vax,infection'*'('*t*','*0.7*','*0.9*')'),
         x="Day") +
    theme_bw() +
    theme(legend.position = "none",
          zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
          zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 
})

plot_alpha1_07_alpha2_09_IIE_vax_death <- lapply(IIE_vax_alpha1_07_alpha2_09, function(df) {
  ggplot(df) +
    geom_line(aes(x = t, y= IIE_death)) +
    labs(y=bquote('IIE'^'Vax,death'*'('*t*','*0.7*','*0.9*')'),
         x="Day") +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal',
          zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
          zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 
})

## (3) collate the plots -------------------------------------------------------
plotlist_eFig5 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_alpha1_07_alpha2_09_IIE_vax_infection$scenario_1,
    i = plot_alpha1_07_alpha2_09_IIE_vax_death$scenario_1 + labs(caption = NULL),
    j = plot_alpha1_07_alpha2_09_IIE_vax_infection$scenario_2,
    k = plot_alpha1_07_alpha2_09_IIE_vax_death$scenario_2 + labs(caption = NULL),
    l = plot_alpha1_07_alpha2_09_IIE_vax_infection$scenario_3,
    m = plot_alpha1_07_alpha2_09_IIE_vax_death$scenario_3 + labs(caption = NULL),
    n = plot_alpha1_07_alpha2_09_IIE_vax_infection$scenario_4,
    o = plot_alpha1_07_alpha2_09_IIE_vax_death$scenario_4 + labs(caption = NULL),
    p = plot_alpha1_07_alpha2_09_IIE_vax_infection$scenario_5,
    q = plot_alpha1_07_alpha2_09_IIE_vax_death$scenario_5 + labs(caption = NULL),
    z = guide_area())

plot_alpha1_07_alpha2_09_IIE_vax_labeled <- wrap_plots(plotlist_eFig5, guides = 'collect', design = layoutplot_eFig3) 

ggsave("3_figures/eFig5.png", plot_alpha1_07_alpha2_09_IIE_vax_labeled, width = 12, height=12, dpi=300, units="in")
