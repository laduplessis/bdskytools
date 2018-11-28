library(bdskytools)

MSE  <- function(x,y) {
  sum((x-y)^2)/length(x)
}

###############################################################################
# 1) Check that density compares to likelihood of trajectory of a single
#    timepoint for likelihood and log-likelihood
#
#   x0    ~ 3
#   mu    ~ 1
#   sigma ~ 1
#   nu    ~ 5
x0 <- 3
mu <- sigma <- 1
nu <- 5

x <- seq(0.5,1.5,length.out=101)
t <- 1
d <- densityOU(x, t, x0, mu, sigma, nu)
l <- c()
for (xx in x) {
    l <- c(l, likOU(c(x0,xx),c(0,t),mu,sigma,nu))
}
print(MSE(d,l))
plot(x,d,type='l',col=pal.dark(cblue))
lines(x,l,lty=2,col=pal.dark(cred))

logd <- densityOU(x, t, x0, mu, sigma, nu, log=TRUE)
logl <- c()
for (xx in x) {
  logl <- c(logl, logLikOU(c(x0,xx),c(0,t),mu,sigma,nu))
}
print(MSE(logd,logl))
plot(x,logd,type='l',col=pal.dark(cblue))
lines(x,logl,lty=2,col=pal.dark(cred))

print(MSE(log(d),logd))
print(MSE(log(l),logl))


# 2) Log-likelihood of a few trajectories, with and without constant term

dt <- 0.1
tn <- 1
t  <- seq(from=0, to=tn, by=dt)

x0    <- 10
mu    <- 5
sigma <- 1
nu    <- 2
x <- simulateOU(x0, t, mu, sigma, nu)
l <- logLikOU(x,t,mu,sigma,nu)
p <- likOU(x,t,mu,sigma,nu)
print(MSE(log(p),l))
print(paste0("{",paste(x,collapse=","),"}  {",paste(t,collapse=","),"}  ",l))

x     <- rep(5,11)
t     <- seq(from=0, to=10, by=1)
x0    <- 5
mu    <- 5
sigma <- 1
nu    <- 2
l     <- logLikOU(x,t,mu,sigma,nu, removeconstant = TRUE)
print(paste0("{",paste(x,collapse=","),"}  {",paste(t,collapse=","),"}  ",l))

x     <- c(10,9,8,7,6,5,5,5,5,6,5)
t     <- seq(from=0, to=10, by=1)
x0    <- 10
mu    <- 5
sigma <- 1
nu    <- 2
l     <- logLikOU(x,t,mu,sigma,nu, removeconstant = TRUE)
print(paste0("{",paste(x,collapse=","),"}  {",paste(t,collapse=","),"}  ",l))

x     <- c(10,9,8,7,6,5,5,5,5,6,5)
t     <- seq(from=0, to=10, by=1)
x0    <- 10
mu    <- 5
sigma <- 1
nu    <- 2
l     <- logLikOU(x,t,mu,sigma,nu, removeconstant = TRUE)
print(paste0("{",paste(x,collapse=","),"}  {",paste(t,collapse=","),"}  ",l))


