#' Simulate BM trajectory
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' 
#' @export
simulateBM <- function(x0, t, k=1) {
  n  <- length(t)
  dt <- diff(t)
  dW <- matrix(0,nrow=k,ncol=n-1)
  W  <- matrix(0,nrow=k,ncol=n)
  for (i in 1:k) {
    dW[i,] <- sqrt(dt)*rnorm(n-1)
    W[i,2:n] <- cumsum(dW[i,])
  }
  
  return(x0+W)
}

#' Calculate quantiles of BM process at different timepoints
#' 
#' @param x0 Initial value of the process 
#' @param t  Vector of time points to evaluate the trajectory at
#' 
#' @export
quantilesBM <- function(p, x0, t) {
  
  E <- x0
  S <- sqrt(t)
  
  result <- matrix(0, nrow=length(p), ncol=length(t))
  for (i in 1:length(p)) {
    result[i, ] <- qnorm(p[i], mean=E, sd=S)
  }
  
  return(result)
}