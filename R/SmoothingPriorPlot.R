
#' Return a function with parameters that can be evaluated using eval function
#' 
#' @export
getPrior <- function(distr, ...) {
  return(substitute(distr(...)))
}


#' Plot the quantiles for OU or BM process as a gradient across the timepoints
#' 
#' @param quantiles The quantiles as determined by quantilesOU or quantilesBM
#' @param t  Vector of time points to plot quantiles at
#' @param col A colour constant, e.g. cblue or cgreen (integer from 1-10)
#' 
#' @export
plotQuantileGradient <- function(t, quantiles, colidx) {
  
  n <- nrow(quantiles)
  col <- pal.dark(colidx, alpha=0.1+0.5*(1/(n/2)))
  
  for (i in 1:floor(n/2)) {
      polygon(c(t, rev(t)), c(quantiles[i,], rev(quantiles[n-i+1,])), col=col, border=NA)
  } 
}

#' Parameters may be a prior or a constant
#' Use getPrior to quote the prior function to be passed with parameters
#' 
#' Should probably rewrite this using do.call()
#' 
#' use par("usr") to get clipping coordinates and then use clip() to set the clipping region
#' 
#' @export
plotOUProcessPriors <- function(x0_prior, t, mu_prior, sigma_prior, nu_prior, ...) {
  
  # Number of traces to simulate
  n <- 10000
  
  # Do simulations
  upper <- lower <- X <-  matrix(0, nrow=n, ncol=length(t))
  for (i in 1:n) {
    x0    <- eval(x0_prior)
    mu    <- eval(mu_prior)
    sigma <- eval(sigma_prior)
    nu    <- eval(nu_prior)
    X[i,] <- simulateOU(x0, t, mu, sigma, nu)
    
    limits    <- quantilesOU(c(0.025, 0.975), x0, t, mu, sigma, nu)
    upper[i,] <- limits[2,]
    lower[i,] <- limits[1,]
  }
  X_hpd <- getMatrixHPD(X)
  limits_hpd <- rbind(getMatrixHPD(lower)[1,], getMatrixHPD(upper)[3,])
  
  
  # Plot traces
  ylims <- paddedrange(limits_hpd)
  plotSkylinePretty(t, X_hpd, type="step", col=pal.dark(cred), fill=pal.dark(cred,0.1), col.axis=pal.dark(cblue),  
                    xaxis=TRUE, yaxis=TRUE, xlab="Time", ylab="OU-Process", xline=2, yline=2, ylims=ylims, ...)
  plotSkylinePretty(t, X, type="lines", traces=max(1000,floor(n/10)), col=pal.dark(cblue,0.1),xaxis=FALSE, yaxis=FALSE, ylims=ylims, new=TRUE, add=TRUE)
  lines(t, limits_hpd[2,], col=pal.light(cred), lty=2)
  lines(t, limits_hpd[1,], col=pal.light(cred), lty=2)
}

