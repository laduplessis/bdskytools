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

#plot(1,type="n",xlim=c(0,tn),ylim=c(0,2),axes=TRUE,xlab=NA,ylab=NA)    
#plotQuantileGradient(t, percentiles, cpurple)
plotOUProcessHPD(x0, t, mu, sigma, nu, ylim=c(0,2), xlab="t",ylab="x")      

for (i in 1:100) {
  x <- simulateOU(x0, t, mu, sigma, nu)
  lines(t, x, col=pal.dark(cpurple,0.25))
}
lines(c(0,tn), rep(mu,2), col=pal.dark(cred), lwd=2)




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
# mu_p    <- getPrior(rnorm, 1, mean=1, sd=0.1)
sigma_p <- getPrior(rlnorm, 1, meanlog=-1, sdlog=1)
# nu_p    <- getPrior(rlnorm, 1, meanlog=-1, sdlog=1)

plotOUProcessHPDEmpirical(x0_p, t, mu, sigma_p, nu, nTraj=1000, plotTraj=TRUE, ylim=c(0,5))
  

####################################################
# Example 4 (rescaling time to get the same process)


#   t2 is scaled up by a factor of 10
#   t3 is same as t1 but with fewer points (necessitates no rescaling)
t1 <- seq(from=0, to=1, by=0.01)
t2 <- seq(from=0, to=10, by=0.1)
t3 <- seq(from=0, to=1, by=0.25)

par(mfrow=c(3,1))

# Theoretical comparison for rescaling time
#   sigma <- sigma/sqrt(factor)
#   nu    <- nu/factor
X_t1 <- plotOUProcessHPD(x0=0.5, t1, mu=1, sigma=0.2, nu=5, ylim=c(0.5,1.2))
X_t2 <- plotOUProcessHPD(x0=0.5, t2, mu=1, sigma=0.2/sqrt(10), nu=5/10, ylim=c(0.5,1.2))
X_t3 <- plotOUProcessHPD(x0=0.5, t3, mu=1, sigma=0.2, nu=5, ylim=c(0.5,1.2))
MSE  <- sum((X_t1-X_t2)^2)/length(X_t1)
print(MSE)

# Empirical comparison for rescaling time
#   sigma <- sigma/sqrt(factor)
#   nu    <- nu/factor
X_t1 <- plotOUProcessHPDEmpirical(x0=0.5, t1, mu=1, sigma=0.2, nu=5, ylim=c(0.5,1.2))
X_t2 <- plotOUProcessHPDEmpirical(x0=0.5, t2, mu=1, sigma=0.2/sqrt(10), nu=5/10, ylim=c(0.5,1.2))
X_t3 <- plotOUProcessHPDEmpirical(x0=0.5, t3, mu=1, sigma=0.2, nu=5, ylim=c(0.5,1.2))
MSE  <- sum((getMatrixHPD(X_t1)-getMatrixHPD(X_t2))^2)/length(X_t1)
print(MSE)


