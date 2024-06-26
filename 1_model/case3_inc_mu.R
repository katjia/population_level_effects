deriv(S_s) <- - beta * S_s * (I_s + I_v) / N 
deriv(S_v) <- - theta * beta * S_v * (I_s + I_v) / N
deriv(I_s) <- beta * S_s * (I_s + I_v) / N - gamma * I_s
deriv(I_v) <- theta * beta * S_v * (I_s + I_v) / N - gamma * I_v 
deriv(R_s) <- (1 - mu) * gamma * I_s 
deriv(R_v) <- (1 - kappa * mu) * gamma * I_v 
deriv(D_s) <- mu * gamma * I_s
deriv(D_v) <- kappa * mu * gamma * I_v

mu <- interpolate(flux_mu_t, flux_mu_y, "constant")
flux_mu_t[] <- user()
flux_mu_y[] <- user()
dim(flux_mu_t) <- user()
dim(flux_mu_y) <- user()
output(mu_t) <- mu 

## Tallying the ever-infected
deriv(Cum_vax) <- 0
deriv(Cum_inf_v) <- theta * beta * S_v * (I_s + I_v) / N
deriv(Cum_inf_s) <- beta * S_s * (I_s + I_v) / N 

## Get N
N <- (S_s + S_v) + (I_s + I_v) + (R_s + R_v)

## Initial states
initial(S_s) <- S_s_ini
initial(S_v) <- S_v_ini

initial(I_s) <- I_s_ini
initial(I_v) <- I_v_ini

initial(R_s) <- R_s_ini
initial(R_v) <- R_v_ini

initial(D_s) <- D_s_ini
initial(D_v) <- D_v_ini

initial(Cum_vax) <- Cum_vax_ini
initial(Cum_inf_v) <- Cum_inf_ini_v
initial(Cum_inf_s) <- Cum_inf_ini_s

## User defined parameters - default in parentheses:
S_s_ini <- user(1000) 
S_v_ini <- user(1000)

I_s_ini <- user(1) 
I_v_ini <- user(0)

R_s_ini <- user(0)
R_v_ini <- user(0)

D_s_ini <- user(0)
D_v_ini <- user(0)

Cum_vax_ini <- user(0)
Cum_inf_ini_v <- user(0)
Cum_inf_ini_s <- user(0)

beta <- user(0.15)  # effective contacts
gamma <- user(0.07) # recovery rate
theta <- user(0.1)  # 1-VE against infection
kappa <- user(0.1)  # 1-VE against death
