lambda <- beta * (I_s + I_v) / N 

# Model parameters
beta <- interpolate(flux_beta_t, flux_beta_y, "constant")
flux_beta_t[] <- user()
flux_beta_y[] <- user()
dim(flux_beta_t) <- user()
dim(flux_beta_y) <- user()
output(beta_t) <- beta

mu <- interpolate(flux_mu_t, flux_mu_y, "constant")
flux_mu_t[] <- user()
flux_mu_y[] <- user()
dim(flux_mu_t) <- user()
dim(flux_mu_y) <- user()
output(mu_t) <- mu 

theta <- interpolate(flux_theta_t, flux_theta_y, "linear")
flux_theta_t[] <- user()
flux_theta_y[] <- user()
dim(flux_theta_t) <- user()
dim(flux_theta_y) <- user()
output(theta_t) <- theta

kappa <- interpolate(flux_kappa_t, flux_kappa_y, "linear")
flux_kappa_t[] <- user()
flux_kappa_y[] <- user()
dim(flux_kappa_t) <- user()
dim(flux_kappa_y) <- user()
output(kappa_t) <- kappa

gamma <- user(0.07) # recovery rate

# Model equations
deriv(S_s) <- - lambda* S_s 
deriv(S_v) <- - theta * lambda * S_v 
deriv(I_s) <- lambda* S_s - gamma * I_s
deriv(I_v) <- theta * lambda* S_v - gamma * I_v 
deriv(R_s) <- (1 - mu) * gamma * I_s
deriv(R_v) <- (1 - kappa * mu) * gamma * I_v
deriv(D_s) <- mu * gamma * I_s
deriv(D_v) <- kappa * mu * gamma * I_v

## Tally the ever-infected
deriv(Cum_vax) <- 0
deriv(Cum_inf_v) <- theta * lambda* S_v
deriv(Cum_inf_s) <- lambda* S_s

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

## User defined initial states
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
