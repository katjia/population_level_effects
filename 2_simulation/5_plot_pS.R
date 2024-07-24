library(RColorBrewer)
# <I> function -----------------------------------------------------------------
get_instant_pS <- function(df_all_counterfactuals){
  df_out <- df_all_counterfactuals %>%
    select(alpha, t, S_s, S_v, N_s, N_v) %>%
    mutate(pS_s = S_s / N_s,
           pS_v = S_v / N_v,
           label = case_when(alpha == 0.9 ~ "90% vaccinated",
                             alpha == 0.7 ~ "70% vaccinated",
                             alpha == 0.0 ~ "0% vaccinated")) 
  return(df_out)
}

# <II> get prop susc -----------------------------------------------------------
df_instant_pS <-
  lapply(df_all_counterfactuals_scenarios_1_to_5, function(df) {
    df %>%
      get_instant_pS()
  })

# <III> plot -------------------------------------------------------------------
plot_instant_pS <- lapply(df_instant_pS, function(df) {
  ggplot(df) +
    geom_line(data = df, aes(x = t, y= pS_s, col="Unvaccinated")) +
    geom_line(data = df, aes(x = t, y= pS_v, col="Vaccinated")) +
    facet_grid(.~label)+
    labs(col="Vaccination status",
         y="Proportion susceptible",
         x="Day") +
    theme_bw() +
    scale_color_manual(values = c("#E78AC3","#66C2A5")) +
    lims(y=c(0,1))
  }) # Expected warning: 731 rows containing missing values due to zeros vaccinated under alpha=0

# collate the plots
layoutplot_eFig6 <- "
ahhhhhhh
bjjjjjjj
clllllll
dnnnnnnn
eppppppp
"

plotlist_eFig6 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    h = plot_instant_pS[["scenario_1"]],
    j = plot_instant_pS[["scenario_2"]],
    l = plot_instant_pS[["scenario_3"]],
    n = plot_instant_pS[["scenario_4"]],
    p = plot_instant_pS[["scenario_5"]])

plot_instant_pS_labeled <- wrap_plots(plotlist_eFig6, guides = 'collect', design = layoutplot_eFig6) 

ggsave("3_figures/eFig6.png", plot_instant_pS_labeled, width = 11, height=9, dpi=300, units="in")

