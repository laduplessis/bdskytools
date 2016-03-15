rm(list = ls())
library(bdskytools)


##########################
# Example 1

set.seed(1)

# Settings for process
nu    <- 5
mu    <- 1
sigma <- 0.7
dt    <- 0.025
tn    <- 1

t  <- seq(from=0, to=tn, by=dt)
x0 <- runif(1, min=0.25, max=1.75)

# 100 percentiles of the process
percentiles <- quantilesOU(seq(from=0, to=1, length.out=101), x0, t, mu, sigma, nu)

# Quantiles for 95% confidence interval and median
limits      <- quantilesOU(c(0.025, 0.5, 0.975), x0, t, mu, sigma, nu)


plot(1,type="n",xlim=c(0,tn),ylim=c(0,2),axes=TRUE,xlab=NA,ylab=NA)    
plotQuantileGradient(t, percentiles, cblue)
for (i in 1:100) {
  x <- simulateOU(x0, t, mu, sigma, nu)
  lines(t, x, col=pal.dark(cgreen,0.25))
}
for (i in 1:3)
  lines(t, limits[i,], col=pal.light(cred), lty=2)
lines(c(0,tn), rep(mu,2), col=pal.light(cred))




##########################
# Example 2

set.seed(10)

n <- 1000

# Settings for process
nu    <- 6
mu    <- 1
sigma <- 0.3
dt    <- 0.025
tn    <- 1

t  <- seq(from=0, to=tn, by=dt)

upper <- lower <- matrix(0, nrow=n, ncol=length(t))

plot(1,type="n",xlim=c(0,tn),ylim=c(0,2),axes=TRUE,xlab=NA,ylab=NA)    
for (i in 1:n) {
  x0 <- runif(1, min=0.25, max=1.75)
  x  <- simulateOU(x0, t, mu, sigma, nu)
  lines(t, x, col=pal.dark(cblue,0.1))
  limits <- quantilesOU(c(0.025, 0.975), x0, t, mu, sigma, nu)
  upper[i,] <- limits[2,]
  lower[i,] <- limits[1,]
}
lines(c(0,tn), rep(mu,2), col=pal.light(cred))
lines(t, getMatrixHPD(upper)[3,], col=pal.light(cred), lty=2)
lines(t, getMatrixHPD(lower)[1,], col=pal.light(cred), lty=2)


##########################
# Example 3

nu    <- 6
mu    <- 1

t    <- seq(from=0, to=1, length.out=40)

x0_p    <- getPrior(rexp, 1, rate=1)
mu_p    <- getPrior(rnorm, 1, mean=1, sd=0.1)
sigma_p <- getPrior(rlnorm, 1, meanlog=-1, sdlog=1)
nu_p    <- getPrior(rlnorm, 1, meanlog=-1, sdlog=1)

plotOUProcessPriors(x0_p, t, mu, sigma_p, nu)
