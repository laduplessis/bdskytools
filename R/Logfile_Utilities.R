###############################################################################
# Utilities for reading BEAST2 logfiles and getting HPDs

#' Wrapper for timing a function
#' (probably not accurate)
#' 
#' @export
timer <- function(fun, ...) {
  t <- Sys.time()
  x <- fun(...)
  print(Sys.time()-t)
  return(x)
}


#' Reverses each column (margin=2) or row (margin=1) of a matrix
revMatrix <- function(data, margin=1) {
  temp <- apply(data, margin, rev)
  if (margin == 1)
    return(t(temp))
  else 
    return(temp)
}


#' Get HPD of a posterior sample
#' Uses Chen and Shao algorithm as implemented in boa package.
#' 
#' @param data The samples from the posterior.
#' @param alpha The confidence level.
#' @return c(lower, median, upper)
#' 
#' @export
getHPD <- function(data, alpha=0.05) {
  hpd <- boa::boa.hpd(data, alpha)
  med <- median(data)
  return(c(hpd[1], med, hpd[2]))
}


#' Get HPD of a matrix of values e.g. a skyline
#' Uses Chen and Shao algorithm as implemented in boa package.
#'  
#' @param data The samples from the posterior. Assumes by default that each row represents a posterior sample and each column a parameter (interval).
#' @param alpha The confidence level.
#' @return 3xn matrix where each column is c(lower, median, upper)
#' 
#' @export
getMatrixHPD <- function(datamat, margin=2, alpha=0.05) {
    return(apply(datamat, margin, getHPD, alpha))
}


#' Extract all matching parameters from the logfile
#' e.g. if par="R0" extract (R0s.1 R0s.2 R0s.3 etc.)
#' 
#' @export
getSkylineSubset <- function(logfile, par) {
  return(lf[grepl(par, names(lf))])
}


#' Read in BEAST2 logfile
#' 
#' @param filename The logfile.
#' @param burnin Discard this percentage of samples.
#' @param maxamples If > 0 stop after reading in this many lines (only for testing).
#' 
#' @export
readLogfile <- function(filename, burnin=0.1, maxsamples=-1) {
  logfile <- read.table(filename, sep="\t", header=TRUE, nrows=maxsamples)
  n <- nrow(logfile)
  return(logfile[floor(burnin*n):n,])
}


#' Interpolates skyline on a time grid (returns a matrix)
#' Assumes that generations are rows and skyline variables are columns (in order)
#' Every generation has the skyline from times[1] to origin[i] (or origin[i] to times[1] if reverse=TRUE)
#' 
#' @param reverse If TRUE assumes that skyline[,n] is the oldest, else if FALSE assume skyline[,1] is the oldest
#' 
#' @export
gridSkyline <- function(skyline, origin, times, reverse=FALSE) {
  
  n <- ncol(skyline)   # skyline variables (shifts+1)
  m <- nrow(skyline)   # generations
  if (length(origin) == 1) origin <- rep(origin,m)
  if (length(origin) != m) stop("Something wrong with dimensions!")
  
  skyline_gridded           <- matrix(0, nrow=m, ncol=length(times))
  colnames(skyline_gridded) <- times
  skyline_matrix            <- as.matrix(skyline)
  
  if (reverse == TRUE)
      getIndices <- function(i) { pmin(n, ceiling(times / origin[i] * n)) }
  else 
      getIndices <- function(i) { pmax(1, n - floor(times / origin[i] * n)) }
  
  for (i in 1:m) {
    #ind    <- pmax(1,n - floor(times / origin[i] * n))
    skyline_gridded[i,] <- skyline_matrix[i,getIndices(i)]
  }
  

  return (skyline_gridded)
}


#' Version using apply, but it is slower for some enigmatic reason
gridSkylineVec <- function(skyline, origin, times, reverse=FALSE) {
  
  n <- ncol(skyline)   # skyline variables (shifts+1)
  m <- nrow(skyline)   # generations
  if (length(origin) != m) stop("Something wrong with dimensions!")
  
  skyline_matrix <- as.matrix(skyline)
  
  getRows <- function(skyline_origin, n, times) {
    ind <- (pmax(1,n - floor(times / skyline_origin[n+1] * n)))
    return(skyline_origin[ind])
  }
  
  skyline_gridded <- apply(cbind(skyline_matrix, origin), 1, getRows, n, times)
  
  return(t(skyline_gridded))
}

