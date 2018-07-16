#' Functions for plotting OU processes
#'



#' Plot expected 95% quantiles (0.025 and 0.975) and median of OU-process, 
#' conditioning on values of x0, mu, sigma, nu
#' 
#' Can also plot a gradient of quantiles (by default plots percentiles - n=101)
#' 
#' 
#' @export
plotOUProcessHPD <- function(x0, t, mu, sigma, nu, 
                             nGrad=101, plotGrad=TRUE, ylim=NULL, ...) {
  
  
  # Get theoretical 95% HPD and median
  limits <- quantilesOU(c(0.025, 0.5, 0.975), x0, t, mu, sigma, nu)
  
  if (is.null(ylim)) {
    ylim <- range(limits)
  }
  plot(1,type='n',xlim=range(t), ylim=ylim, ...)
  
  # Plot gradient of quantiles
  if (plotGrad == TRUE) { 
    percentiles <- quantilesOU(seq(from=0, to=1, length.out=nGrad), x0, t, mu, sigma, nu)
    plotQuantileGradient(t, percentiles, cblue)
  }
  
  # Plot 95% HPD and median
  lines(t, limits[1,], col=pal.dark(cred), lty=2, lwd=2)
  lines(t, limits[2,], col=pal.dark(cred), lty=2, lwd=2)
  lines(t, limits[3,], col=pal.dark(cred), lty=2, lwd=2)
  
  return(limits)
}


#' Plot empirical 95% quantiles (0.025 and 0.975) and median of OU-process, 
#' from simulating nTraj replicates, which can also be plotted
#' 
#' Instead of specifying constant parameter values, prior functions can also
#' be specified, in which case they will be sampled for simulating trajectories
#' 
#' @export
plotOUProcessHPDEmpirical <- function(x0_prior, t, mu_prior, sigma_prior, nu_prior, 
                                      nTraj=1000, plotTraj=TRUE, ylim=c(0,5), ...) {
  
  # Setup
  X <-  matrix(0, nrow=nTraj, ncol=length(t))
  plot(1,type="n",xlim=range(t), ylim=ylim,axes=TRUE,xlab=NA,ylab=NA, ...)    
  
  # Do simulations and plot trajectories
  for (i in 1:nTraj) {
    x0    <- eval(x0_prior)
    mu    <- eval(mu_prior)
    sigma <- eval(sigma_prior)
    nu    <- eval(nu_prior)
    
    # Simulate trajectory
    X[i,] <- simulateOU(x0, t, mu, sigma, nu)
    
    if (plotTraj) lines(t, X[i,], col=pal.dark(cblue,0.1), lwd=0.5)
  }
  
  # Empirical 95% HPD and median
  X_hpd <- getMatrixHPD(X)
  lines(t, X_hpd[1,], col=pal.dark(cred), lty=1, lwd=2)
  lines(t, X_hpd[2,], col=pal.dark(cred), lty=1, lwd=2)
  lines(t, X_hpd[3,], col=pal.dark(cred), lty=1, lwd=2)
  
  return(X)
  
}