#' Functions for simulating Ornstein-Uhlenbeck processes and calculating densities and likelihoods
#' 


#' Simulate OU trajectory using the Euler-Maruyama method (inaccurate)
#' 
#' Function assumes equidistant time intervals
#'
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
#' @export
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
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
#' @export
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
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
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


#' Expected value of OU process at times t starting at x0 at time 0
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
#' @export
expectedValOU <- function(x0, t, mu, nu) {
  return(mu + (x0-mu)*exp(-nu*t))
}

#' Standard deviation of OU process at times t starting at x0 at time 0
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
#' @export
standardDevOU <- function(x0, t, mu, sigma, nu) {
  term <- 1 - exp(-2*nu*t)
  return(sqrt((sigma^2/(2*nu))*term))
}


#' Calculate the density of OU process for a given (x,t) pair
#' 
#' @param x  Vector of trajectory points to calculate density of
#' @param t  Vector of time points associated with x
#' @param x0 Initial value of the process 
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' @param log Return log(p)
#' 
#' @export
densityOU <- function(x,t,x0,mu,sigma,nu,log=FALSE) {
  
  E <- expectedValOU(x0, t, mu, nu)
  S <- standardDevOU(x0, t, mu, sigma, nu)
  
  d <- dnorm(x, mean=E, sd = S, log=log)
  
  return(d)
}

#' Calculate quantiles of OU process at different timepoints
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
#' @export
quantilesOU <- function(p, x0, t, mu, sigma, nu) {
  
  E <- expectedValOU(x0, t, mu, nu)
  S <- standardDevOU(x0, t, mu, sigma, nu)
  
  result <- matrix(0, nrow=length(p), ncol=length(t))
  for (i in 1:length(p)) {
    result[i, ] <- qnorm(p[i], mean=E, sd=S)
  }
  
  return(result)
}



#' Calculates the likelihood of observing a trajectory under a given
#' OU-process
#' 
#' Likelihood calculation follows from memorylessness of a Markov chain, 
#' so it's just a product of Gaussians
#' 
#' @param x = (x0, x1, ..., xn)
#' @param t = (t0, t2, ..., tn)
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' 
#' @export
likOU <- function(x,t,mu,sigma,nu) {
  
  if (nu < 0 || sigma < 0)    
      return(0)
  
  n <- length(x)
  L <- 1
  
  for (i in 2:n) {
    dt <- t[i]-t[i-1]
    
    meanterm <- mu + (x[i-1] - mu)*exp(-nu*dt)
    varterm  <- (sigma^2)*(1 - exp(-2*nu*dt))/(2*nu)
    
    sqterm   <- (x[i]-meanterm)^2
    
    L <- L * exp(-sqterm/(2*varterm))/sqrt(2*pi*varterm) 
  }
  return(L)
}

#' Calculates the log-likelihood of observing a trajectory under a given
#' OU-process
#' 
#' Likelihood calculation follows from memorylessness of a Markov chain, 
#' so it's just a sum of Gaussians
#' 
#' @param x = (x0, x1, ..., xn)
#' @param t = (t0, t2, ..., tn)
#' @param mu Mean
#' @param sigma Standard deviation
#' @param nu Decay rate
#' @param removeconstant Drop the constant term
#' 
#' @export
logLikOU <- function(x,t,mu,sigma,nu,removeconstant=FALSE) {
  
  if (nu < 0 || sigma < 0)    
    return(-Inf)
  
  n    <- length(x)
  logP <- - 0.5*(n-1)*log((sigma^2)/(2*nu))
  
  if (removeconstant == FALSE) 
    logP <- logP -0.5*(n-1)*log(2*pi)
  
  for (i in 2:n) {
    dt <- t[i]-t[i-1]
    
    relterm  <- 1 - exp(-2*nu*dt)
    diffterm <- x[i] - mu - (x[i-1] - mu)*exp(-nu*dt)
    
    term <- (diffterm^2)/relterm
    logP <- logP - 0.5*log(relterm) - (nu/(sigma^2))*term
  }
  
  return(logP)
}


