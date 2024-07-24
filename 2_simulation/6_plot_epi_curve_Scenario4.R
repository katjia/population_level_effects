# <I> plot epidemic curves under Scenario 4 ------------------------------------
df_scenario_4 <- df_all_counterfactuals_scenarios_1_to_5[["scenario_4"]] %>%
  filter(alpha %in% c(0.7,0.9)) %>%
  mutate(eff_susc = theta_t*S_v + S_s) 

sir_curves_scenario_4 <- df_scenario_4 %>%
  select(t, "1. Susceptible" = S, "2. Infectious"= I, "3. Recovered" = R, alpha) %>%
  pivot_longer(-c(t, alpha), values_to="Count", names_to = "compartment") 

epi_curves_sceario_4 <- ggplot() +
  geom_line(data=sir_curves_scenario_4, aes(y=Count, x=t, col=compartment)) +
  geom_line(data=df_scenario_4, aes(y=eff_susc, x=t), lty="dashed") +
  facet_grid(.~alpha) +
  labs(x="Day",
       col="Compartment")+
  theme_bw() 

ggsave("3_figures/eFig7.png", epi_curves_sceario_4, width = 6.5, height=3.5, dpi=300, units="in")

