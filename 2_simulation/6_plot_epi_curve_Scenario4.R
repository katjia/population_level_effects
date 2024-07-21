# Plot epidemic curve under scenario 4
df_scenario4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0.7,0.9)) %>%
  mutate(eff_susc = theta_t*S_v + S_s) 

sir_curve_scenario4 <- df_scenario4 %>%
  select(t, "1. Susceptible" = S, "2. Infectious"= I, "3. Recovered" = R, alpha) %>%
  pivot_longer(-c(t, alpha), values_to="Count", names_to = "compartment") 

epi_curve_sceario4 <- ggplot() +
  geom_line(data=sir_curve_scenario4, aes(y=Count, x=t, col=compartment)) +
  geom_line(data=df_scenario4, aes(y=eff_susc, x=t), lty="dashed") +
  facet_grid(.~alpha) +
  labs(x="Day",
       col="Compartment")+
  theme_bw() 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig8.png", epi_curve_sceario4, width = 6.5, height=3.5, dpi=300, units="in")

