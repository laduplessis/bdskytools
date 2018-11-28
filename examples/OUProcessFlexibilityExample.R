library(bdskytools)
rm(list = ls())

x0   <- c(1,1.5,2,5)
xlim <- c(0,5)

#####################################################################
# Example 1: With fixed parameter values (get theoretical quantiles)

plotOUProcessFlexibility(xlim, dt=0.025, x0, mu=1, sigma=1, nu=1) 
plotOUProcessFlexibility(xlim, dt=0.1,   x0, mu=1, sigma=1, nu=1) 
plotOUProcessFlexibility(xlim, dt=10,    x0, mu=1, sigma=1, nu=1) 

plotOUProcessFlexibility(xlim, dt=0.025, x0, mu=1, sigma=1, nu=10) 
plotOUProcessFlexibility(xlim, dt=0.1,   x0, mu=1, sigma=1, nu=10) 
plotOUProcessFlexibility(xlim, dt=10,    x0, mu=1, sigma=1, nu=10) 

plotOUProcessFlexibility(xlim, dt=0.025, x0, mu=1, sigma=5, nu=10) 
plotOUProcessFlexibility(xlim, dt=0.1,   x0, mu=1, sigma=5, nu=10) 
plotOUProcessFlexibility(xlim, dt=10,    x0, mu=1, sigma=5, nu=10) 


#####################################################################
# Example 2: Sample parameters from priors 
# (get empirical HPDs from simulated trajectories)

mu_p    <- getPrior(rlnorm, 1, meanlog=0, sdlog=1.25)
sigma_p <- getPrior(rlnorm, 1, meanlog=0, sdlog=0.25)
nu_p    <- getPrior(rexp, 1, rate=1)

plotOUProcessFlexibilityEmpirical(xlim, dt=0.025, x0, mu_p, sigma_p, nu_p) 
plotOUProcessFlexibilityEmpirical(xlim, dt=0.1, x0, mu_p, sigma_p, nu_p) 
plotOUProcessFlexibilityEmpirical(xlim, dt=0.5, x0, mu_p, sigma_p, nu_p) 

