library(ggforce)

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
    summarise(IIE_vax_infection = last(pInf_v) - first(pInf_v),
              IIE_vax_death = last(pDeath_v) - first(pDeath_v))
}

# <I> Plot IIE^Unvax(t, alpha0=0, alpha1=0.7) for Claim 1a ---------------------------
## (1) Case 1: Time-invariant parameters ----------------------------------------
### (1.1) get IIE 
df_alpha0_0_alpha1_07_IIE_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_IIE_unvax()

### (1.2) plot 
plot_alpha0_0_alpha1_07_IIE_infection_case1 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case1, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha0_0_alpha1_07_IIE_death_case1 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case1, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

## (2) Case 2: Increasing effective contacts (beta) -----------------------------
### (2.1) get IIE
df_alpha0_0_alpha1_07_IIE_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_IIE_unvax()

### (2.2) plot 
plot_alpha0_0_alpha1_07_IIE_infection_case2 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case2, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha0_0_alpha1_07_IIE_death_case2 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case2, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25))

## (3) Case 3: Increasing IFR (mu) ----------------------------------------------
### (3.1) get IIE
df_alpha0_0_alpha1_07_IIE_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_IIE_unvax()

## (3.2) plot 
plot_alpha0_0_alpha1_07_IIE_infection_case3 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case3, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha0_0_alpha1_07_IIE_death_case3 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case3, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

## (4) Case 4: Waning VEs ------------------------------------------------------
### (4.1) get IE
df_alpha0_0_alpha1_07_IIE_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_IIE_unvax()

### (4.2) plot Case 4: waning VEs 
plot_alpha0_0_alpha1_07_IIE_infection_case4 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case4, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha0_0_alpha1_07_IIE_death_case4 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case4, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25))

## (5) Case 5: Cases 2 and 4 -----------------------------------
### (5.1) get IIE
df_alpha0_0_alpha1_07_IIE_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0, 0.7)) %>%
  get_IIE_unvax()

### (5.3) plot
plot_alpha0_0_alpha1_07_IIE_infection_case5 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case5, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha0_0_alpha1_07_IIE_death_case5 <- ggplot() +
  geom_line(data = df_alpha0_0_alpha1_07_IIE_case5, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0*','*0.7*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 5", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_eFig4 <- "
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

plotlist_eFig4 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_alpha0_0_alpha1_07_IIE_infection_case1,
    i = plot_alpha0_0_alpha1_07_IIE_death_case1 + labs(caption = NULL),
    j = plot_alpha0_0_alpha1_07_IIE_infection_case2,
    k = plot_alpha0_0_alpha1_07_IIE_death_case2 + labs(caption = NULL),
    l = plot_alpha0_0_alpha1_07_IIE_infection_case3,
    m = plot_alpha0_0_alpha1_07_IIE_death_case3 + labs(caption = NULL),
    n = plot_alpha0_0_alpha1_07_IIE_infection_case4,
    o = plot_alpha0_0_alpha1_07_IIE_death_case4 + labs(caption = NULL),
    p = plot_alpha0_0_alpha1_07_IIE_infection_case5,
    q = plot_alpha0_0_alpha1_07_IIE_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_alpha0_0_alpha1_07_IIE_labeled <- wrap_plots(plotlist_eFig4, guides = 'collect', design = layoutplot_eFig4) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig4.png", plot_alpha0_0_alpha1_07_IIE_labeled, width = 12, height=12, dpi=300, units="in")


# <II> Plot IIE^Unax(t, alpha1=0.7, alpha2=0.9) for Claim 1b ------------------
## (1) Case 1: Time-invariant parameters ---------------------------------------
### (1.1) get IIE 
df_alpha1_07_alpha2_09_IIE_unvax_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_unvax()

### (1.2) plot 
plot_alpha1_07_alpha2_09_IIE_unvax_infection_case1 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case1, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_unvax_death_case1 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case1, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

## (2) Case 2: Increasing effective contacts (beta) -----------------------------
### (2.1) get IIE
df_alpha1_07_alpha2_09_IIE_unvax_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_unvax()

### (2.2) plot 
plot_alpha1_07_alpha2_09_IIE_unvax_infection_case2 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case2, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_unvax_death_case2 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case2, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25))

## (3) Case 3: Increasing IFR (mu) ----------------------------------------------
### (3.1) get IIE
df_alpha1_07_alpha2_09_IIE_unvax_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_unvax()

## (3.2) plot 
plot_alpha1_07_alpha2_09_IIE_unvax_infection_case3 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case3, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_unvax_death_case3 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case3, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

## (4) Case 4: Waning VEs ------------------------------------------------------
### (4.1) get IE
df_alpha1_07_alpha2_09_IIE_unvax_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_unvax()

### (4.2) plot Case 4: waning VEs 
plot_alpha1_07_alpha2_09_IIE_unvax_infection_case4 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case4, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_unvax_death_case4 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case4, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25))

## (5) Case 5: Cases 2 and 4 -----------------------------------
### (5.1) get IIE
df_alpha1_07_alpha2_09_IIE_unvax_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_unvax()

### (5.3) plot
plot_alpha1_07_alpha2_09_IIE_unvax_infection_case5 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case5, aes(x = t, y= IIE_infection)) +
  labs(y=bquote('IIE'^'Unvax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_unvax_death_case5 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_unvax_case5, aes(x = t, y= IIE_death)) +
  labs(y=bquote('IIE'^'Unvax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 5", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_eFig5 <- "
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

plotlist_eFig5 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_alpha1_07_alpha2_09_IIE_unvax_infection_case1,
    i = plot_alpha1_07_alpha2_09_IIE_unvax_death_case1 + labs(caption = NULL),
    j = plot_alpha1_07_alpha2_09_IIE_unvax_infection_case2,
    k = plot_alpha1_07_alpha2_09_IIE_unvax_death_case2 + labs(caption = NULL),
    l = plot_alpha1_07_alpha2_09_IIE_unvax_infection_case3,
    m = plot_alpha1_07_alpha2_09_IIE_unvax_death_case3 + labs(caption = NULL),
    n = plot_alpha1_07_alpha2_09_IIE_unvax_infection_case4,
    o = plot_alpha1_07_alpha2_09_IIE_unvax_death_case4 + labs(caption = NULL),
    p = plot_alpha1_07_alpha2_09_IIE_unvax_infection_case5,
    q = plot_alpha1_07_alpha2_09_IIE_unvax_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_alpha1_07_alpha2_09_IIE_unvax_labeled <- wrap_plots(plotlist_eFig5, guides = 'collect', design = layoutplot_eFig5) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig5.png", plot_alpha1_07_alpha2_09_IIE_unvax_labeled, width = 12, height=12, dpi=300, units="in")



# <III> Plot IIE^Vax(t, alpha1=0.7, alpha2=0.9) for Claim 1b ------------------
## (1) Case 1: Time-invariant parameters ---------------------------------------
### (1.1) get IIE 
df_alpha1_07_alpha2_09_IIE_vax_case1 <- df_all_counterfactuals_case1 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_vax()

### (1.2) plot 
plot_alpha1_07_alpha2_09_IIE_vax_infection_case1 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case1, aes(x = t, y= IIE_vax_infection)) +
  labs(y=bquote('IIE'^'Vax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_vax_death_case1 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case1, aes(x = t, y= IIE_vax_death)) +
  labs(y=bquote('IIE'^'Vax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

## (2) Case 2: Increasing effective contacts (beta) -----------------------------
### (2.1) get IIE
df_alpha1_07_alpha2_09_IIE_vax_case2 <- df_all_counterfactuals_case2 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_vax()

### (2.2) plot 
plot_alpha1_07_alpha2_09_IIE_vax_infection_case2 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case2, aes(x = t, y= IIE_vax_infection)) +
  labs(y=bquote('IIE'^'Vax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_vax_death_case2 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case2, aes(x = t, y= IIE_vax_death)) +
  labs(y=bquote('IIE'^'Vax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25))

## (3) Case 3: Increasing IFR (mu) ----------------------------------------------
### (3.1) get IIE
df_alpha1_07_alpha2_09_IIE_vax_case3 <- df_all_counterfactuals_case3 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_vax()

## (3.2) plot 
plot_alpha1_07_alpha2_09_IIE_vax_infection_case3 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case3, aes(x = t, y= IIE_vax_infection)) +
  labs(y=bquote('IIE'^'Vax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_vax_death_case3 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case3, aes(x = t, y= IIE_vax_death)) +
  labs(y=bquote('IIE'^'Vax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

## (4) Case 4: Waning VEs ------------------------------------------------------
### (4.1) get IE
df_alpha1_07_alpha2_09_IIE_vax_case4 <- df_all_counterfactuals_case4 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_vax()

### (4.2) plot Case 4: waning VEs 
plot_alpha1_07_alpha2_09_IIE_vax_infection_case4 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case4, aes(x = t, y= IIE_vax_infection)) +
  labs(y=bquote('IIE'^'Vax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_vax_death_case4 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case4, aes(x = t, y= IIE_vax_death)) +
  labs(y=bquote('IIE'^'Vax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25))

## (5) Case 5: Cases 2 and 4 -----------------------------------
### (5.1) get IIE
df_alpha1_07_alpha2_09_IIE_vax_case5 <- df_all_counterfactuals_case5 %>%
  filter(alpha %in% c(0.7, 0.9)) %>%
  get_IIE_vax()

### (5.3) plot
plot_alpha1_07_alpha2_09_IIE_vax_infection_case5 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case5, aes(x = t, y= IIE_vax_infection)) +
  labs(y=bquote('IIE'^'Vax,infection'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = "none",
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

plot_alpha1_07_alpha2_09_IIE_vax_death_case5 <- ggplot() +
  geom_line(data = df_alpha1_07_alpha2_09_IIE_vax_case5, aes(x = t, y= IIE_vax_death)) +
  labs(y=bquote('IIE'^'Vax,death'*'('*t*','*0.7*','*0.9*')'),
       x="Day") +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        zoom.x = element_rect(fill = "white", color = "black", linewidth = 0.25), 
        zoom.y = element_rect(fill = "white", color = "black", linewidth = 0.25)) 

# collate the plots
row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 5", angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death") + theme_void() 

layoutplot_eFig6 <- "
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

plotlist_eFig6 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    e = row5,
    f = col1,
    g = col2,
    h = plot_alpha1_07_alpha2_09_IIE_vax_infection_case1,
    i = plot_alpha1_07_alpha2_09_IIE_vax_death_case1 + labs(caption = NULL),
    j = plot_alpha1_07_alpha2_09_IIE_vax_infection_case2,
    k = plot_alpha1_07_alpha2_09_IIE_vax_death_case2 + labs(caption = NULL),
    l = plot_alpha1_07_alpha2_09_IIE_vax_infection_case3,
    m = plot_alpha1_07_alpha2_09_IIE_vax_death_case3 + labs(caption = NULL),
    n = plot_alpha1_07_alpha2_09_IIE_vax_infection_case4,
    o = plot_alpha1_07_alpha2_09_IIE_vax_death_case4 + labs(caption = NULL),
    p = plot_alpha1_07_alpha2_09_IIE_vax_infection_case5,
    q = plot_alpha1_07_alpha2_09_IIE_vax_death_case5 + labs(caption = NULL),
    z = guide_area())

plot_alpha1_07_alpha2_09_IIE_vax_labeled <- wrap_plots(plotlist_eFig6, guides = 'collect', design = layoutplot_eFig6) 

ggsave("~/Documents/GitHub/population_level_effects/3_figures/eFig6.png", plot_alpha1_07_alpha2_09_IIE_vax_labeled, width = 12, height=12, dpi=300, units="in")
