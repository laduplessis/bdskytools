###############################################################################
# Functions for plotting prior distributions (also OU priors)


#' Plot a parametric distribution where the distribution is already available in R
#' 
#' For the distribution "xxx" the functions "qxxx()" and "dxxx()" need to be defined and available
#' 
#' @param priorfun Family of distributions e.g. norm or beta
#' @param prior_args List of arguments for priorfun e.g. list(mean=1,sd=0.5)
#' @param meaninrealspace Only applies to lognormal distribution
#' 
#' @export
plotPrior <- function(priorfun="norm", prior_args=list(), col=pal.dark(cblue), fill=pal.dark(cblue, 0.25), ylab="Density", xlab="x", 
                      plotquantile=TRUE, grid=TRUE, positive=TRUE, invert=FALSE, meaninrealspace=FALSE, scaling=1, 
                      add=FALSE, xlims=NULL, ylims=NULL, ...) {

  
  scaleX <- function(x) {
      if (invert == TRUE)
          return(scaling/x)
      else
          return(scaling*x)
  }
  
  
  if (meaninrealspace == TRUE && priorfun == "lnorm") {
      prior_args[[1]] <- log(prior_args[[1]])-0.5*prior_args[[2]]^2
  }
    
  
  # Get 95% HPD
  q <- do.call(sprintf("q%s",priorfun), c(list(p=c(0.025, 0.5, 0.975)), prior_args))
  names(q) <- c(0.025, 0.5, 0.975)
  print(paste("95% quantiles =",paste(q,collapse = ", ")))
  
  # Get limits
  if (is.null(xlims)) {
      xlims <- paddedrange(q)
      if (positive == TRUE)
          xlims[1] <- max(0, xlims[1])
  }
  
  x <- seq(xlims[1], xlims[2], length.out=1000)
  dfun  <- do.call(sprintf("d%s",priorfun), c(list(x=x), prior_args))
  
  # Get quantiles as indices on x
  q_idx <- sapply(q, function(i) { max(which(x <= i)) })
  
  
  if (is.null(ylims)) {
      ylims <- paddedrange(dfun)
      ylims[1] <- 0
  }
  
  # Deal with values at Infinity
  dfun[which(dfun == Inf)] <- ylims[2]*10
    
  # Plot the prior
  # Do not open a new device (add to the current plot)
  if (add == FALSE) 
      plot(1,type='n',xlim=sort(scaleX(xlims)), ylim=ylims, las=1, axes=TRUE, ylab=ylab, xlab=xlab, ...)

  if (grid==TRUE) { 
      for (y in axTicks(2)) 
        abline(h=y, lty=3, lwd=0.5)
  }
  
  # Draw distribution
  
  # Left tail
  polygon(scaleX(c(xlims[1], x[1:q_idx[1]], x[q_idx[1]])), c(ylims[1], dfun[1:q_idx[1]], ylims[1]), 
          col=pal.dark(cred,0.25), border=NA)
  lines(scaleX(x[1:q_idx[1]]), dfun[1:q_idx[1]], col=pal.dark(cred), lwd=2)
  
  # 95%
  polygon(scaleX(c(x[q_idx[1]], x[q_idx[1]:q_idx[3]], x[q_idx[3]])), c(ylims[1], dfun[q_idx[1]:q_idx[3]], ylims[1]), 
          col=fill, border=NA)
  lines(scaleX(x[q_idx[1]:q_idx[3]]), dfun[q_idx[1]:q_idx[3]], col=col, lwd=2)
  
  # Right tail
  polygon(scaleX(c(x[q_idx[3]], x[q_idx[3]:length(x)], xlims[2])), c(ylims[1], dfun[q_idx[3]:length(dfun)], ylims[1]), 
          col=pal.dark(cred,0.25), border=NA)
  lines(scaleX(x[q_idx[3]:length(x)]), dfun[q_idx[3]:length(x)], col=pal.dark(cred), lwd=2)
  
  # Draw median
  #lines(rep(x[q_idx[2]],2), c(ylims[1],dfun[q_idx[2]]), col=col, lwd=2)
  abline(v=scaleX(q[2]), lty=2, lwd=2 ,col=pal.dark(cred))
}


# Paratially inspired by blog post from Jordan Douglas (https://www.beast2.org/2017/05/05/using-r-to-estimate-probability-distribution-parameters-using-quantiles-or-hpd-intervals/)
optimPrior <- function(priorfun="norm", prior_probs=c(0.025, 0.5, 0.975), prior_quantiles=c(-1.959964, 0, 1.959964), prior_init=c(1,1)) {
  
  pfun     <- function(x,params) do.call(sprintf("p%s",priorfun), c(list(x), params))
  errorfun <- function(params) sum((sapply(prior_quantiles,pfun,params)-prior_probs)^2)
  result   <- optim(prior_init, errorfun)
  
  return(result)
}


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

#' This function does not really plot the true quantiles - quantiles of quantiles across sims
#' (function should be removed)
#' 
#' Parameters may be a prior or a constant
#' Use getPrior to quote the prior function to be passed with parameters
#' 
#' Should probably rewrite this using do.call()
#' 
#' @export
plotOUProcessPriors <- function(x0_prior, t, mu_prior, sigma_prior, nu_prior, ...) {
  
  warning("Be careful! This function does not do what it claims to...")
  
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


#' This function does not really plot the true quantiles - quantiles of quantiles across sims
#' (function should be removed)
#' 
#' Parameters may be a prior or a constant
#' 
#' Use getPrior to quote the prior function to be passed with parameters
#' 
#' 
#' @export
plotBMProcessPriors <- function(x0_prior, t, ...) {
  
  warning("Be careful! This function does not do what it claims to...")
  
  
  # Number of traces to simulate
  n <- 10000
  
  # Do simulations
  upper <- lower <- X <-  matrix(0, nrow=n, ncol=length(t))
  for (i in 1:n) {
    x0    <- eval(x0_prior)
    X[i,] <- simulateBM(x0, t, mu)
    
    limits    <- quantilesBM(c(0.025, 0.975), x0, t)
    upper[i,] <- limits[2,]
    lower[i,] <- limits[1,]
  }
  X_hpd <- getMatrixHPD(X)
  limits_hpd <- rbind(getMatrixHPD(lower)[1,], getMatrixHPD(upper)[3,])
  
  
  # Plot traces
  ylims <- paddedrange(limits_hpd)
  plotSkylinePretty(t, X_hpd, type="step", col=pal.dark(cred), fill=pal.dark(cred,0.1), col.axis=pal.dark(cblue),  
                    xaxis=TRUE, yaxis=TRUE, xlab="Time", ylab="BM-Process", xline=2, yline=2, ylims=ylims, ...)
  plotSkylinePretty(t, X, type="lines", traces=max(1000,floor(n/10)), col=pal.dark(cblue,0.1),xaxis=FALSE, yaxis=FALSE, ylims=ylims, new=TRUE, add=TRUE)
  lines(t, limits_hpd[2,], col=pal.light(cred), lty=2)
  lines(t, limits_hpd[1,], col=pal.light(cred), lty=2)
}
