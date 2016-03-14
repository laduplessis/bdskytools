

#' Simulate OU trajectory using the Euler-Maruyama method (inaccurate)
#' 
simulateDiscreteOU <- function(x0, t, mu, sigma, nu) {
  
  dt   <- t[2]-t[1]
  x    <- matrix(data = 0, nrow=1, ncol=length(t))
  x[1] <- x0
  for (i in 1:(length(t)-1)) {
    x[i+1] <- x[i] + nu*(mu-x[i])*dt + sigma*sqrt(dt)*rnorm(1)
  }
  return(x)
}


#' Simulate OU trajectory using the analytical solution obtained
#' using a scaled time-transformed Wiener process
simulateOU <- function(x0, t, mu, sigma, nu) {
  
  n  <- length(t)
  dt <- diff(t)
  dW <- sqrt(diff(exp(2*nu*t)-1))*rnorm(n-1)
  W  <- matrix(0,nrow=1,ncol=n)
  ex <- exp(-nu*t)
  
  W[2:n] <- cumsum(dW)
  X      <- x0*ex + mu*(1-ex) + sigma*ex*W/sqrt(2*nu)
  return(X)
}


#' Simulate OU trajectory using the analytical solution obtained
#' using a scaled time-transformed Wiener process
#' 
#' Slow unvectorized version
simulateOUslow <- function(x0, t, mu, sigma, nu) {
  
  dW   <- rnorm(length(t)-1)
  
  W <- matrix(0, nrow=1, ncol=length(t))
  for (i in 1:(length(t)-1)) {
    W[i+1] <- W[i] + sqrt(exp(2*nu*t[i+1])-exp(2*nu*t[i]))*dW[i]
  }
  ex <- exp(-nu*t)
  x <- x0*ex + mu*(1-ex) + sigma*ex*W/sqrt(2*nu)
  
  return(x) 
}


#' Expected value of OU process at times t starting 
#' at x0 at time 0
expectedValOU <- function(x0, t, mu, nu) {
  return(mu + (x0-mu)*exp(-nu*t))
}

#' Standard deviation of OU process at times t starting 
#' at x0 at time 0
standardDevOU <- function(x0, t, mu, sigma, nu) {
  term <- 1 - exp(-2*nu*t)
  return(sqrt((sigma^2/(2*nu))*term))
}

#' Calculate quantiles of OU process at different timepoints
#' (assuming x0 is the only information)
quantilesOU <- function(p, x0, t, mu, sigma, nu) {
  
  E <- expectedValOU(x0, t, mu, nu)
  S <- standardDevOU(x0, t, mu, sigma, nu)
  
  result <- matrix(0, nrow=length(p), ncol=length(t))
  for (i in 1:length(p)) {
    result[i, ] <- qnorm(p[i], mean=E, sd=S)
  }
  
  return(result)
}


plotQuantileGradient <- function(quantiles, col) {
  
  n <- nrow(quantiles)
  col <- pal.dark(col, alpha=0.1+0.5*(1/(n/2))) 
  
  for (i in 1:floor(n/2)) {
    polygon(c(t, rev(t)), c(quantiles[i,], rev(quantiles[n-i+1,])), col=col, border=NA)
  }
  
}
