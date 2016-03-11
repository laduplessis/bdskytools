
#' Reverses each column (margin=2) or row (margin=1) of a matrix
revMatrix <- function(data, margin=1) {
  temp <- apply(data, margin, rev)
  if (margin == 1)
    return(t(temp))
  else 
    return(temp)
}

#' Interpolates skyline on a time grid (returns a matrix)
#' Assumes that generations are rows and skyline variables are columns (in order)
#' Every generation has the skyline from times[1] to origin[i] (or origin[i] to times[1] if reverse=TRUE)
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
  
  for (i in 1:m) {
    ind    <- pmax(1,n - floor(times / origin[i] * n))
    skyline_gridded[i,] <- skyline_matrix[i,ind]
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


#' Get range and pad by some percentage of the range on either side
#' Good for plot limits
paddedrange <- function(y, pad=0.1) {
  yrange <- c(min(y), max(y))
  return( c(yrange[1]-diff(yrange)*pad*0.5, yrange[2]+diff(yrange)*pad*0.5) )
}

plotSkylineSmooth <- function(hpdmatrix, times, col=pal.dark(cblack), fill=pal.dark(cgray, 0.25), ...) {

  ylims <- paddedrange(hpdmatrix)
  xlims <- range(times)

  plot(1, type='n', ylim=ylims, xlim=xlims, ...)
  polygon(c(times, rev(times)), c(hpdmatrix[1,], rev(hpdmatrix[3,])), col=fill, border=NA)
  lines(times, hpdmatrix[2,], col=col)
}


plotSkylineStepped <- function(hpdmatrix, times,  col=pal.dark(cblack), fill=pal.dark(cgray, 0.25), ...) {
  
  ylims <- paddedrange(hpdmatrix)
  xlims <- range(times)
  
  plot(1, type='n', ylim=ylims, xlim=xlims, ...)
  for (i in 2:ncol(hpdmatrix)) {
      rect(times[i-1], hpdmatrix[1,i-1], times[i], hpdmatrix[3,i-1], col=fill, border=NA)
  }
  lines(times, hpdmatrix[2,], col=col, type='S')
}

plotSkylineTraces <- function(skyline_mat, times, traces=1000, col=pal.dark(cgray, 0.1), type='s', ...) {
  
  ylims <- paddedrange(skyline_mat)
  xlims <- range(times)

  if (traces > 0 && traces < nrow(skyline_mat)) 
    ind <- sample(nrow(skyline_mat), size=traces, replace=FALSE)
  else 
    ind <- 1:nrow(skyline_mat)
      
  plot(1, type='n', ylim=ylims, xlim=xlims, ...)
  for (i in 1:length(ind)) {
    lines(times, skyline_mat[ind[i],], col=col, type=type)
  }
}


#' new: create a new set of axes
#' add: add to the current plot (do not create a new plotting device)
#' 
#' assume each column is a different time point (ncol(skylinemat) == length(times) must be fulfilled)
#' 
#' @export
plotSkyline <- function(skyline_mat, times, type="", col=pal.dark(cblack), fill=pal.dark(cgray, 0.25), new=TRUE, add=FALSE, ...) {
  
  # Do not open a new device (add to the current plot)
  if (add == TRUE) par(new=TRUE)
  
  # Create a new set of axes that fits this skyline
  if (new == TRUE) {
      ylims <- paddedrange(skyline_mat)
      xlims <- range(times)
      
      plot(1, type='n', ylim=ylims, xlim=xlims, ...)
  }
  
  polygon(c(times, rev(times)), c(skyline_mat[1,], rev(skyline_mat[3,])), col=fill, border=NA)
  lines(times, skyline_mat[2,], col=col)
  
  par(new=FALSE)
}


