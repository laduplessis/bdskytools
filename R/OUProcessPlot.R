#######################################
# Functions for plotting OU processes #
#######################################
# Louis du Plessis, 2018              #
#######################################


#' Plot expected 95% quantiles (0.025 and 0.975) and median of OU-process, 
#' conditioning on values of x0, mu, sigma, nu
#' 
#' Can also plot a gradient of quantiles (by default plots percentiles - n=101)
#' 
#' TODO: Name is misleading, quantiles are NOT HPDs.
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


#' Plot empirical 95% HPD interval and median of OU-process, 
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
  plot(1,type="n",xlim=range(t), ylim=ylim,axes=TRUE, ...)    
  
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




#' Plot the density of an OU-process with parameters x0, mu, sigma, nu after time-interval
#' dt, at the range of values given by xrange.
#' 
#' Density is only plotted between the 0.025 and 0.975 quantiles.
#' 
#' @export
plotOUDensity <- function(xrange, dt, x0, mu, sigma, nu, col=pal.dark(cblue), col.alpha=0.25, 
                          new=TRUE, plotMedian=FALSE, plotX0=TRUE, plotGrid=TRUE, ...) { 
  
  OUdensity <- densityOU(xrange,dt,x0,mu,sigma,nu)
  OUlimits  <- quantilesOU(c(0.025, 0.5, 0.975), x0, dt, mu, sigma, nu)
  CIlimits  <- findInterval(OUlimits, xrange)
  
  if (new) {
    plot(c(xrange[1],xrange), c(0,OUdensity), type='n', ...)
    if (plotGrid) {
        grid(col='black')
    }
  }
  
  lines(xrange[CIlimits[1]:CIlimits[3]], OUdensity[CIlimits[1]:CIlimits[3]], lwd=2, col=col, xpd=TRUE)
  polygon(c(xrange[max(1,CIlimits[1])],xrange[CIlimits[1]:CIlimits[3]],xrange[CIlimits[3]]), c(0,OUdensity[CIlimits[1]:CIlimits[3]],0), col=scales::alpha(col,col.alpha), border=NA)
  
  if (plotX0) {
    points(x0,0,pch=16,cex=2,col=col,xpd=TRUE)
  }
  
  if (plotMedian) {
    abline(v=OUlimits[2], col=col, lwd=2, lty=2)
  }
  
  return(c(x0=x0, shift=(OUlimits[2]-x0), width=unname(OUlimits[3]-OUlimits[1]), centering=(OUlimits[2]-OUlimits[1])/(OUlimits[3]-OUlimits[1])))
}


#' Plot the density of ntraj OU-processes simulated to time dt, using priors
#' x0_prior, mu_prior, sigma_prior and nu_prior on the OU-process parameters
#' (fixed values are also supported).
#' 
#' Density is only plotted in the 95% HPD interval, at the points in xrange.
#'
#' @export
plotOUDensityEmpirical <- function(xrange, dt, x0_prior, mu_prior, sigma_prior, nu_prior, nTraj=1000, nSteps=100, col=pal.dark(cblue), col.alpha=0.25, 
                                   bw='sj', new=TRUE, plotMedian=FALSE, plotX0=TRUE, plotGrid=TRUE, ...) { 
  
  
  # Setup
  t <- seq(0,dt,length.out=nSteps)
  X <- matrix(0, nrow=nTraj, ncol=length(t))
  
  # Do simulations and plot trajectories
  x0list <- c()
  for (i in 1:nTraj) {
    x0    <- eval(x0_prior)
    mu    <- eval(mu_prior)
    sigma <- eval(sigma_prior)
    nu    <- eval(nu_prior)
    
    # Simulate trajectory
    X[i,] <- simulateOU(x0, t, mu, sigma, nu)
    
    x0list <- c(x0list, x0)
  }
  # Empirical 95% HPD and median
  OUlimits  <- getHPD(X[,ncol(X)])
  OUdensity <- density(X[,ncol(X)], bw=bw, from=OUlimits[1], to=OUlimits[3])
  
  if (new) {
    plot(c(xrange[1],OUdensity$x,xrange[length(xrange)]), c(0,OUdensity$y,0), type='n', ...)
    if (plotGrid) {
        grid(col='black')
    }
  }
  
  lines(OUdensity$x, OUdensity$y, lwd=2, col=col)
  polygon(c(OUdensity$x[1], OUdensity$x, OUdensity$x[length(OUdensity$x)]), c(0, OUdensity$y, 0), col=scales::alpha(col,col.alpha), border=NA)
  
  if (plotX0) {
    points(x0,0,pch=16,cex=2,col=col,xpd=TRUE)
  }
  
  if (plotMedian) {
    abline(v=OUlimits[2], col=col, lwd=2, lty=2)
  }
  
  
  return(c(x0=median(x0list), shift=(OUlimits[2]-median(x0list)), width=unname(OUlimits[3]-OUlimits[1]), centering=(OUlimits[2]-OUlimits[1])/(OUlimits[3]-OUlimits[1])))
}


#' Plot a summary of the flexibility of an OU-process with parameters mu, sigma and nu, starting 
#' from different x0 values, over a time period of dt
#' 
#' Answers the questions: 
#' - How big is the 95% CI after one time interval? 
#  - How much can the mean revert after one time interval?
#' 
#' @export
plotOUProcessFlexibility <- function(xlim, dt, x0, mu, sigma, nu, nSteps=100) { 
  
  t    <- seq(0,dt, by=dt/nSteps)
  x    <- seq(xlim[1],xlim[2], length.out=1000)
  n    <- length(x0)
  cols <- RColorBrewer::brewer.pal(max(3,n),"Set2")
  
  layout(matrix(c(1:n,rep(n+1,n)),byrow=TRUE, nrow=2))
  for (i in 1:n) {
    plotOUProcessHPD(x0[i], t, mu, sigma, nu, ylim=range(x), xlab="t",ylab="x",main=paste("x0 =",x0[i]))
  }
  
  result <- c()
  for (i in 1:n) {
    result <- rbind(result, plotOUDensity(x, dt, x0[i], mu, sigma, nu, col=cols[i], new=ifelse(i==1,TRUE,FALSE),xaxs='i', yaxs='i', xlab='x', ylab="", bty='n'))
  } 
  abline(v=mu,col='red',lty=3, lwd=2)
  legend('topright',col=cols, pch=16, bty='n', legend=paste("x0 =",x0), cex=1.25)
  
  return(result)
}




#' Plot a summary of the flexibility of an OU-process with parameters sampled from priors 
#' mu_prior, sigma_prior and nu_prior, starting from different x0 values, over a time period 
#' of dt 
#' 
#' Answers the questions: 
#' - How big is the 95% CI after one time interval? 
#  - How much can the mean revert after one time interval?
#' 
#' @export
plotOUProcessFlexibilityEmpirical <- function(xlim, dt, x0, mu_p, sigma_p, nu_p, nSteps=100) { 
  
  t    <- seq(0,dt, by=dt/nSteps)
  x    <- seq(xlim[1],xlim[2], length.out=1000)
  n    <- length(x0)
  cols <- RColorBrewer::brewer.pal(max(3,n),"Set2")
  
  layout(matrix(c(1:n,rep(n+1,n)),byrow=TRUE, nrow=2))
  for (i in 1:n) {
    plotOUProcessHPDEmpirical(x0[i], t, mu_p, sigma_p, nu_p, nTraj=1000, plotTraj=TRUE, ylim=range(x), xlab="t",ylab="x", main=paste("x0 =",x0[i]))
  }
  
  result <- c()
  for (i in 1:n) {
    result <- rbind(result, plotOUDensityEmpirical(x, dt, x0[i], mu_p, sigma_p, nu_p, col=cols[i], new=ifelse(i==1,TRUE,FALSE),xaxs='i', yaxs='i', xlab='x', ylab="", bty='n'))
  } 
  
  mu <- c()
  for (i in 1:5000) {
    mu <- c(mu, eval(mu_p))
  }
  #print(length(mu))
  abline(v=median(mu),col='red',lty=3, lwd=2)
  legend('topright',col=cols, pch=16, bty='n', legend=paste("x0 =",x0), cex=1.25)
  
  return(result)
}



###############################################################################

testFlexibility <- function() {
  
  x0   <- c(1,1.5,2,5)
  dt   <- 0.5
  xlim <- c(0,8)
  
  nu    <- 1
  mu    <- 1
  sigma <- 1
  plotOUProcessFlexibility(xlim, dt, x0, mu, sigma, nu) 
  
  mu_p    <- getPrior(rlnorm, 1, meanlog=0, sdlog=1.25)
  sigma_p <- getPrior(rlnorm, 1, meanlog=0, sdlog=0.25)
  nu_p    <- getPrior(rexp, 1, rate=1)
  plotOUProcessFlexibilityEmpirical(xlim, dt, x0, mu_p, sigma_p, nu_p) 
}

