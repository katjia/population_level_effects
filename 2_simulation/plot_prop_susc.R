library(RColorBrewer)
library(patchwork)

# (1) Case 1: Time-invariant parameters ----------------------------------------
## get prop susceptible
### a function to get instantaneous proportion susceptible
get_instant_pS <- function(dataframe){
  dataframe %>%
    filter(VE_infection %in% c(0, 0.5, 0.9, 1)) %>%
    select(VE_infection, alpha, t, S_s, S_v, N_s, N_v, N) %>%
    mutate(pS_s = S_s / N_s,
           pS_v = S_v / N_v,
           label = case_when(alpha == 0.7 ~ "70% vaccinated",
                             alpha == 0.0 ~ "0% vaccinated")) 
}

df_instant_pS_case1 <- get_instant_pS(df_all_counterfactuals_case1_time_invariant_par)

plot_instant_pS_case1 <- ggplot() +
  geom_line(data = df_instant_pS_case1, aes(x = t, y= pS_s, col="Unvaccinated", lty = as.factor(VE_infection))) +
  geom_line(data = df_instant_pS_case1, aes(x = t, y= pS_v, col="Vaccinated", lty = as.factor(VE_infection))) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%"))+
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# (2) Case 2: Increasing effective contacts (beta) -----------------------------
df_instant_pS_case2 <- get_instant_pS(df_all_counterfactuals_case2_inc_beta)

plot_instant_pS_case2 <- ggplot() +
  geom_line(data = df_instant_pS_case2, aes(x = t, y= pS_s, col="Unvaccinated", lty = as.factor(VE_infection))) +
  geom_line(data = df_instant_pS_case2, aes(x = t, y= pS_v, col="Vaccinated", lty = as.factor(VE_infection))) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%"))+
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# (3) Case 3: Increasing IFR (mu) ----------------------------------------------
df_instant_pS_case3 <- get_instant_pS(df_all_counterfactuals_case3_inc_mu)

plot_instant_pS_case3 <- ggplot() +
  geom_line(data = df_instant_pS_case3, aes(x = t, y= pS_s, col="Unvaccinated", lty = as.factor(VE_infection))) +
  geom_line(data = df_instant_pS_case3, aes(x = t, y= pS_v, col="Vaccinated", lty = as.factor(VE_infection))) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%"))+
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# (4) Case 4: Waning VEs -------------------------------------------------------
df_instant_pS_case4 <- get_instant_pS(df_all_counterfactuals_case4_waning) 

plot_instant_pS_case4 <- ggplot() +
  geom_line(data = df_instant_pS_case4, aes(x = t, y= pS_s, col="Unvaccinated", lty = as.factor(VE_infection))) +
  geom_line(data = df_instant_pS_case4, aes(x = t, y= pS_v, col="Vaccinated", lty = as.factor(VE_infection))) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%"))+
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

## (5) Case 5: Cases 2 and 4 ---------------------------------------------------
df_instant_pS_case5 <- get_instant_pS(df_all_counterfactuals_case5_inc_beta_waning) 

plot_instant_pS_case5 <- ggplot() +
  geom_line(data = df_instant_pS_case5, aes(x = t, y= pS_s, col="Unvaccinated", lty = as.factor(VE_infection))) +
  geom_line(data = df_instant_pS_case5, aes(x = t, y= pS_v, col="Vaccinated", lty = as.factor(VE_infection))) +
  facet_grid(.~label)+
  labs(col="Vaccination status",
       y="Proportion susceptible",
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%"))+
  theme_bw() +
  scale_color_manual(values = c("#E78AC3","#66C2A5")) +
  lims(y=c(0,1))

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 5", angle = 90) + theme_void() 

layoutplot3 <- "
ahhhhhhh
bjjjjjjj
clllllll
dnnnnnnn
eppppppp
"

plotlist3 <-
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

plot_instant_pS_labeled <- wrap_plots(plotlist3, guides = 'collect', design = layoutplot3) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig3_traj_pS.png", plot_instant_pS_labeled, width = 10, height=9, dpi=300, units="in")

