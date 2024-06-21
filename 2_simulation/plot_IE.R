library(ggforce)
# (1) Case 1: Time-invariant parameters ---------------------------------------
## (1.1) compute pY and pD
df_pY_pD_case1 <- df_all_counterfactuals_case1_time_invariant_par %>%
  filter(alpha %in% c(0, 0.7)) %>%
  filter(VE_infection %in% c(0, 0.5, 0.9, 1)) %>%
  select(-par, VE_infection, alpha, t, Cum_inf_s, D_s) %>%
  mutate(
    pInf_s = case_when(alpha == 0.7 ~ Cum_inf_s / (2e4 * 0.3),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
    pDeath_s = case_when(alpha == 0.7 ~ D_s / (2e4 * 0.3),
                         alpha == 0.0 ~ D_s / 2e4)
  )
## (1.2) compute IE
df_IE_death_case1 <- df_pY_pD_case1 %>%
  select(VE_infection, alpha, t, pDeath_s) %>%
  group_by(VE_infection, t) %>%
  pivot_wider(names_from = alpha, values_from = pDeath_s) %>%
  mutate(IE_death = `0` - `0.7`)

df_IE_infection_case1 <- df_pY_pD_case1 %>%
  select(VE_infection, alpha, t, pInf_s) %>%
  group_by(VE_infection, t) %>%
  pivot_wider(names_from = alpha, values_from = pInf_s) %>%
  mutate(IE_infection = `0` - `0.7`)

## (1.3) plot 
plot_IE_infection_case1 <- ggplot() +
  geom_line(data = df_IE_infection_case1, aes(x = t, y= IE_infection, lty = as.factor(VE_infection))) +
  labs(tag = "A(i)",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.01,0.05), zoom.size = .75) +
  theme_bw() +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"), 
                        breaks=c(0, 0.5, 0.9, 1), 
                          labels = c("0%, 90%", 
                                     "50%, 90%", 
                                     "90%, 90%", 
                                     "100%, 90%")) +
  theme(legend.position = "none")

plot_IE_death_case1 <- ggplot() +
  geom_line(data = df_IE_death_case1, aes(x = t, y= IE_death, lty = as.factor(VE_infection))) +
  labs(tag = "A(ii)",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*0.7*')'),
       x="Day",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case1_time_invariant_par.R","\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_IE.R")) +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                     "50%, 90%", 
                                     "90%, 90%", 
                                     "100%, 90%")) +
  facet_zoom(ylim = c(-0.001,0.005), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

# (2) Case 2: Increasing effective contacts (beta) -----------------------------
## (2.1) compute pY and pD
df_pY_pD_case2 <- df_all_counterfactuals_case2_inc_beta %>%
  filter(alpha %in% c(0, 0.7)) %>%
  filter(VE_infection %in% c(0, 0.5, 0.9, 1)) %>%
  select(-par, VE_infection, alpha, t, Cum_inf_s, D_s) %>%
  mutate(pInf_s = case_when(alpha == 0.7 ~ Cum_inf_s / (2e4 * 0.3),
                          alpha == 0.0 ~ Cum_inf_s / 2e4),
         pDeath_s = case_when(alpha == 0.7 ~ D_s / (2e4 * 0.3),
                          alpha == 0.0 ~ D_s / 2e4))
## (2.2) compute IE
df_IE_death_case2 <- df_pY_pD_case2 %>%
  select(VE_infection, alpha, t, pDeath_s) %>%
  group_by(VE_infection, t) %>%
  pivot_wider(names_from = alpha, values_from = pDeath_s) %>%
  mutate(IE_death = `0` - `0.7`)

df_IE_infection_case2 <- df_pY_pD_case2 %>%
  select(VE_infection, alpha, t, pInf_s) %>%
  group_by(VE_infection, t) %>%
  pivot_wider(names_from = alpha, values_from = pInf_s) %>%
  mutate(IE_infection = `0` - `0.7`)

## (2.3) plot 
plot_IE_infection_case2 <- ggplot() +
  geom_line(data = df_IE_infection_case2, aes(x = t, y= IE_infection, lty = as.factor(VE_infection))) +
  labs(tag="B(i)",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.3,0.3), zoom.size = .75) +
  theme_bw() +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%")) +
  theme(legend.position = "none") 


plot_IE_death_case2 <- ggplot() +
  geom_line(data = df_IE_death_case2, aes(x = t, y= IE_death, lty = as.factor(VE_infection))) +
  labs(tag="B(ii)",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*0.7*')'),
       x="Day",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case2_inc_beta.R","\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_IE.R")) +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%")) +
  facet_zoom(ylim = c(-0.003,0.003), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

# (3) Case 3: Increasing IFR (mu) ----------------------------------------------
## (3.1) compute pY and pD
df_pY_pD_case3 <- df_all_counterfactuals_case3_inc_mu %>%
  filter(alpha %in% c(0, 0.7)) %>%
  filter(VE_infection %in% c(0, 0.5, 0.9, 1)) %>%
  select(-par, VE_infection, alpha, t, Cum_inf_s, D_s) %>%
  mutate(pInf_s = case_when(alpha == 0.7 ~ Cum_inf_s / (2e4 * 0.3),
                          alpha == 0.0 ~ Cum_inf_s / 2e4),
         pDeath_s = case_when(alpha == 0.7 ~ D_s / (2e4 * 0.3),
                          alpha == 0.0 ~ D_s / 2e4))
## (3.2) compute IE
df_IE_death_case3 <- df_pY_pD_case3 %>%
  select(VE_infection, alpha, t, pDeath_s) %>%
  group_by(VE_infection, t) %>%
  pivot_wider(names_from = alpha, values_from = pDeath_s) %>%
  mutate(IE_death = `0` - `0.7`)

df_IE_infection_case3 <- df_pY_pD_case3 %>%
  select(VE_infection, alpha, t, pInf_s) %>%
  group_by(VE_infection, t) %>%
  pivot_wider(names_from = alpha, values_from = pInf_s) %>%
  mutate(IE_infection = `0` - `0.7`)

## (3.3) plot 
plot_IE_infection_case3 <- ggplot() +
  geom_line(data = df_IE_infection_case3, aes(x = t, y= IE_infection, lty = as.factor(VE_infection))) +
  labs(tag="C(i)",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  facet_zoom(ylim = c(-0.3,0.3), zoom.size = .75) +
  theme_bw() +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                                    labels = c("0%, 90%", 
                                               "50%, 90%", 
                                               "90%, 90%", 
                                               "100%, 90%")) +
  theme(legend.position = "none")

plot_IE_death_case3 <- ggplot() +
  geom_line(data = df_IE_death_case3, aes(x = t, y= IE_death, lty = as.factor(VE_infection))) +
  labs(tag="C(ii)",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*0.7*')'),
       x="Day",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case3_inc_mu.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_IE.R")) +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%")) +
  facet_zoom(ylim = c(-0.003,0.003), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

plot_IE_infection_case3 + plot_IE_death_case3

# (4) Case 4: Waning VEs ------------------------------------------------------
## (4.1) compute pY and pD
df_pY_pD_case4 <- df_all_counterfactuals_case4_waning %>%
  filter(alpha %in% c(0, 0.7)) %>%
  filter(VE_infection_t0 %in% c(0, 0.5, 0.9, 1)) %>%
  select(-par, VE_infection_t0, alpha, t, Cum_inf_s, D_s) %>%
  mutate(
    pInf_s = case_when(alpha == 0.7 ~ Cum_inf_s / (2e4 * 0.3),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
    pDeath_s = case_when(alpha == 0.7 ~ D_s / (2e4 * 0.3),
                         alpha == 0.0 ~ D_s / 2e4)
  )
## (4.2) compute IE
df_IE_death_case4 <- df_pY_pD_case4 %>%
  select(VE_infection_t0, alpha, t, pDeath_s) %>%
  group_by(VE_infection_t0, t) %>%
  pivot_wider(names_from = alpha, values_from = pDeath_s) %>%
  mutate(IE_death = `0` - `0.7`)

df_IE_infection_case4 <- df_pY_pD_case4 %>%
  select(VE_infection_t0, alpha, t, pInf_s) %>%
  group_by(VE_infection_t0, t) %>%
  pivot_wider(names_from = alpha, values_from = pInf_s) %>%
  mutate(IE_infection = `0` - `0.7`)

## (4.3) plot Case 4: waning VEs 
plot_IE_infection_case4 <- ggplot() +
  geom_line(data = df_IE_infection_case4, aes(x = t, y= IE_infection, lty = as.factor(VE_infection_t0))) +
  labs(tag="D(i)",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%")) +
  facet_zoom(ylim = c(-0.01,0.05), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = "none")

plot_IE_death_case4 <- ggplot() +
  geom_line(data = df_IE_death_case4, aes(x = t, y= IE_death, lty = as.factor(VE_infection_t0))) +
  labs(tag="D(ii)",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*0.7*')'),
       x="Day",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case4_waning_VEs.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_IE.R")) +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%")) +
  facet_zoom(ylim = c(-0.001,0.005), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

plot_IE_infection_case4 + plot_IE_death_case4

# (5) Case 5: Cases 2 and 4 -----------------------------------
## (5.1) compute pY and pD
df_pY_pD_case5 <- df_all_counterfactuals_case5_inc_beta_waning %>%
  filter(alpha %in% c(0, 0.7)) %>%
  filter(VE_infection_t0 %in% c(0, 0.5, 0.9, 1)) %>%
  select(-par, VE_infection_t0, alpha, t, Cum_inf_s, D_s) %>%
  mutate(
    pInf_s = case_when(alpha == 0.7 ~ Cum_inf_s / (2e4 * 0.3),
                       alpha == 0.0 ~ Cum_inf_s / 2e4),
    pDeath_s = case_when(alpha == 0.7 ~ D_s / (2e4 * 0.3),
                         alpha == 0.0 ~ D_s / 2e4)
  )
## compute IE
df_IE_death_case5 <- df_pY_pD_case5 %>%
  select(VE_infection_t0, alpha, t, pDeath_s) %>%
  group_by(VE_infection_t0, t) %>%
  pivot_wider(names_from = alpha, values_from = pDeath_s) %>%
  mutate(IE_death = `0` - `0.7`)

df_IE_infection_case5 <- df_pY_pD_case5 %>%
  select(VE_infection_t0, alpha, t, pInf_s) %>%
  group_by(VE_infection_t0, t) %>%
  pivot_wider(names_from = alpha, values_from = pInf_s) %>%
  mutate(IE_infection = `0` - `0.7`)

## plot
plot_IE_infection_case5 <- ggplot() +
  geom_line(data = df_IE_infection_case5, aes(x = t, y= IE_infection, lty = as.factor(VE_infection_t0))) +
  labs(tag="E(i)",
       y=bquote('IIE'^'Inf'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%")) +
  facet_zoom(ylim = c(-0.3,0.3), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = "none")

plot_IE_death_case5 <- ggplot() +
  geom_line(data = df_IE_death_case5, aes(x = t, y= IE_death, lty = as.factor(VE_infection_t0))) +
  labs(tag="E(ii)",
       y=bquote('IIE'^'Death'*'('*t*','*0*','*0.7*')'),
       x="Day",
       caption = paste0("~/Documents/GitHub/population_level_effects/1_model/case5_inc_beta_waning_VEs.R", "\n", 
                        "~/Documents/GitHub/population_level_effects/2_simulation/simulation.R", "\n",
                        "~/Documents/GitHub/population_level_effects/2_simulation/plot_IE.R")) +
  scale_linetype_manual(name="VE infection, VE death", 
                        values=c("solid", "22", "42", "73"),  
                        labels = c("0%, 90%", 
                                   "50%, 90%", 
                                   "90%, 90%", 
                                   "100%, 90%"))+
  facet_zoom(ylim = c(-0.003,0.003), zoom.size = .75) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Case 5", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot <- "
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

plotlist <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_IE_infection_case1,
    i = plot_IE_death_case1 + labs(caption = NULL),
    j = plot_IE_infection_case2,
    k = plot_IE_death_case2 + labs(caption = NULL),
    l = plot_IE_infection_case3,
    m = plot_IE_death_case3 + labs(caption = NULL),
    n = plot_IE_infection_case4,
    o = plot_IE_death_case4 + labs(caption = NULL),
    p = plot_IE_infection_case5,
    q = plot_IE_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_IIE_labeled <- wrap_plots(plotlist, guides = 'collect', design = layoutplot) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/Fig1_IIEs.png", plot_IIE_labeled, width = 12, height=12, dpi=300, units="in")
