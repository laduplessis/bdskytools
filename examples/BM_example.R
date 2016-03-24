#rm(list = ls())
library(bdskytools)


##########################
# Example 1

set.seed(1)

# Settings for process
dt <- 0.025
tn <- 1
t  <- seq(from=0, to=tn, by=dt)
x0 <- 1

# 100 percentiles of the process
percentiles <- quantilesBM(seq(from=0, to=1, length.out=101), x0, t)

# Quantiles for 95% confidence interval and median
limits      <- quantilesBM(c(0.025, 0.5, 0.975), x0, t)

# Simulate traces
x <- simulateBM(x0, t,100)

plot(1,type="n",xlim=c(0,tn),ylim=x0+c(-3,3),axes=TRUE,xlab=NA,ylab=NA)    
plotQuantileGradient(t, percentiles, cblue)
for (i in 1:100) {
  lines(t, x[i,], col=pal.dark(corange,0.25))
}
for (i in 1:3)
  lines(t, limits[i,], col=pal.light(cred), lty=2)
lines(c(0,tn), rep(x0,2), col=pal.light(cred))



##########################
# Example 2

set.seed(10)

n <- 1000

# Settings for process
dt    <- 0.025
tn    <- 1

t  <- seq(from=0, to=tn, by=dt)

upper <- lower <- matrix(0, nrow=n, ncol=length(t))

plot(1,type="n",xlim=c(0,tn),ylim=c(-2,4),axes=TRUE,xlab=NA,ylab=NA)    
for (i in 1:n) {
  x0 <- runif(1, min=0.25, max=1.75)
  x  <- simulateBM(x0, t)
  lines(t, x, col=pal.dark(cblue,0.1))
  limits <- quantilesBM(c(0.025, 0.975), x0, t)
  upper[i,] <- limits[2,]
  lower[i,] <- limits[1,]
}
lines(c(0,tn), rep(mu,2), col=pal.light(cred))
lines(t, getMatrixHPD(upper)[3,], col=pal.light(cred), lty=2)
lines(t, getMatrixHPD(lower)[1,], col=pal.light(cred), lty=2)


##########################
# Example 3

t    <- seq(from=0, to=1, length.out=40)

x0_p    <- getPrior(rexp, 1, rate=1)

NewFig("~/Documents/Projects/BEAST2/bdskytools/examples/BMTest.pdf", width=7)
plotBMProcessPriors(x0_p, t)
dev.off()

