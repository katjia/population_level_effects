library(RColorBrewer)
library(patchwork)

# (1) Case 1: Time-invariant parameters ----------------------------------------
### a function to get instantaneous proportion susceptible
get_instant_pS <- function(df_all_counterfactuals){
  df_out <- df_all_counterfactuals %>%
    select(VE_infection, alpha, t, S_s, S_v, N_s, N_v, N) %>%
    mutate(pS_s = S_s / N_s,
           pS_v = S_v / N_v,
           label = case_when(alpha == 0.9 ~ "90% vaccinated",
                             alpha == 0.7 ~ "70% vaccinated",
                             alpha == 0.0 ~ "0% vaccinated")) 
  return(df_out)
}

df_instant_pS_case1 <- get_instant_pS(df_all_counterfactuals_case1)

plot_instant_pS_case1 <- ggplot() +
  geom_line(data = df_instant_pS_case1, aes(x = t, y= pS_s, col="Unvaccinated")) +
  geom_line(data = df_instant_pS_case1, aes(x = t, y= pS_v, col="Vaccinated")) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# (2) Case 2: Increasing effective contacts (beta) -----------------------------
df_instant_pS_case2 <- get_instant_pS(df_all_counterfactuals_case2)

plot_instant_pS_case2 <- ggplot() +
  geom_line(data = df_instant_pS_case2, aes(x = t, y= pS_s, col="Unvaccinated")) +
  geom_line(data = df_instant_pS_case2, aes(x = t, y= pS_v, col="Vaccinated")) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# (3) Case 3: Increasing IFR (mu) ----------------------------------------------
df_instant_pS_case3 <- get_instant_pS(df_all_counterfactuals_case3)

plot_instant_pS_case3 <- ggplot() +
  geom_line(data = df_instant_pS_case3, aes(x = t, y= pS_s, col="Unvaccinated")) +
  geom_line(data = df_instant_pS_case3, aes(x = t, y= pS_v, col="Vaccinated")) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# (4) Case 4: Waning VEs -------------------------------------------------------
df_instant_pS_case4 <- get_instant_pS(df_all_counterfactuals_case4) 

plot_instant_pS_case4 <- ggplot() +
  geom_line(data = df_instant_pS_case4, aes(x = t, y= pS_s, col="Unvaccinated")) +
  geom_line(data = df_instant_pS_case4, aes(x = t, y= pS_v, col="Vaccinated")) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

## (5) Case 5: Cases 2 and 4 ---------------------------------------------------
df_instant_pS_case5 <- get_instant_pS(df_all_counterfactuals_case5) 

plot_instant_pS_case5 <- ggplot() +
  geom_line(data = df_instant_pS_case5, aes(x = t, y= pS_s, col="Unvaccinated")) +
  geom_line(data = df_instant_pS_case5, aes(x = t, y= pS_v, col="Vaccinated")) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 5", angle = 90) + theme_void() 

layoutplot_eFig7 <- "
ahhhhhhh
bjjjjjjj
clllllll
dnnnnnnn
eppppppp
"

plotlist_eFig7 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    h = plot_instant_pS_case1,
    j = plot_instant_pS_case2,
    l = plot_instant_pS_case3,
    n = plot_instant_pS_case4,
    p = plot_instant_pS_case5)

plot_instant_pS_labeled <- wrap_plots(plotlist_eFig7, guides = 'collect', design = layoutplot_eFig7) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig7.png", plot_instant_pS_labeled, width = 11, height=9, dpi=300, units="in")

