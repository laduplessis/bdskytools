library(bdskytools)

###############################################################################
# Make results reproducible
set.seed(10)

# Time intervals
dt <- 0.005
tn <- 1
t  <- seq(from=0, to=tn, by=dt)

par(mfrow=c(3,2))

#######################################
# Process 1 (fixed values to do a sanity check with theoretical quantiles)
#   x0    ~ 3
#   mu    ~ 1
#   sigma ~ 1
#   nu    ~ 5
x0 <- 3
mu <- sigma <- 1
nu <- 5

plotOUProcessHPDEmpirical(x0, t, mu, sigma, nu, ylim=c(0,4), main="Fixed vals sanity check")
limits <- quantilesOU(c(0.025, 0.5, 0.975), x0, t, mu, sigma, nu)
lines(t, limits[1,], col=pal.dark(cgreen), lty=2, lwd=2)
lines(t, limits[2,], col=pal.dark(cgreen), lty=2, lwd=2)
lines(t, limits[3,], col=pal.dark(cgreen), lty=2, lwd=2)
legend("topright", lty=c(1,2), col=pal.dark(c(cred,cgreen)), legend=c("Empirical HPD","Theoretical HPD"))


#######################################
# Process 2
#   x0    ~ unif(0.5,1.5)
#   mu    ~ lognorm(0,0.5)
#   sigma ~ norm(0.5,0.1)
#   nu    ~ gamma(1,5)
x0_p    <- getPrior(runif,  1, min=1.5, max=2.5)     
mu_p    <- getPrior(rlnorm, 1, meanlog=0, sdlog=0.5)
sigma_p <- getPrior(rnorm,  1, mean=0.5, sd=0.1)
nu_p    <- getPrior(rgamma, 1, shape=1, scale=5)

plotOUProcessHPDEmpirical(x0_p, t, mu_p, sigma_p, nu_p, ylim=c(0,4), main="Priors test")



#######################################
# Process 3 (fixed mean, high variance, low mean reversion)
#   x0    ~ exp(1,1)
#   mu    ~ 1
#   sigma ~ norm(1,0.5)
#   nu    ~ gamma(1,5)
x0_p    <- getPrior(rexp, 1, rate=1)
mu_p    <- 1
sigma_p <- getPrior(rnorm,  1, mean=1, sd=0.5)
nu_p    <- getPrior(rgamma, 1, shape=1, scale=5)

plotOUProcessHPDEmpirical(x0_p, t, mu_p, sigma_p, nu_p, ylim=c(0,4), main="Fixed mean, high var, low mean reversion")


#######################################
# Process 4 (fixed mean, low variance, low  mean reversion)
#   x0    ~ exp(1,1)
#   mu    ~ 1
#   sigma ~ norm(0.1,0.1)
#   nu    ~ gamma(1,5)
x0_p    <- getPrior(rexp, 1, rate=1)
mu_p    <- 1
sigma_p <- getPrior(rnorm,  1, mean=0.1, sd=0.1)
nu_p    <- getPrior(rgamma, 1, shape=1, scale=5)

plotOUProcessHPDEmpirical(x0_p, t, mu_p, sigma_p, nu_p, ylim=c(0,4), main="Fixed mean, low var, low mean reversion")



#######################################
# Process 5 (fixed mean, high variance, high mean reversion)
#   x0    ~ exp(1,1)
#   mu    ~ 1
#   sigma ~ norm(1,0.5)
#   nu    ~ gamma(1,50)
x0_p    <- getPrior(rexp, 1, rate=1)
mu_p    <- 1
sigma_p <- getPrior(rnorm,  1, mean=1, sd=0.5)
nu_p    <- getPrior(rgamma, 1, shape=1, scale=50)

plotOUProcessHPDEmpirical(x0_p, t, mu_p, sigma_p, nu_p, ylim=c(0,4), main="Fixed mean, high var, high mean reversion")


#######################################
# Process 6 (fixed mean, low variance, high mean reversion)
#   x0    ~ exp(1,1)
#   mu    ~ 1
#   sigma ~ norm(0.1,0.1)
#   nu    ~ gamma(1,50)
x0_p    <- getPrior(rexp, 1, rate=1)
mu_p    <- 1
sigma_p <- getPrior(rnorm,  1, mean=0.1, sd=0.1)
nu_p    <- getPrior(rgamma, 1, shape=1, scale=50)

plotOUProcessHPDEmpirical(x0_p, t, mu_p, sigma_p, nu_p, ylim=c(0,4), main="Fixed mean, low var, high mean reversion")


