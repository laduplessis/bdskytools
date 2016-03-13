


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


#' Plot a Skyline.
#' 
#' Assume each column is a different time point (ncol(skylinemat) == length(times) must be fulfilled)
#' 
#' @param type Type of skyline to plot.
#' 
#'             "smooth" : If skyline_mat contains the HPDs plot a smooth line and polygon (type='l')
#'             
#'             "step"   : If skyline_mat contains the HPDs plot a stepped line and polygon (type='S')
#'             
#'             "lines" : Plot every row in the matrix as a smooth line (type='l')
#'             
#'             "steplines" : Plot every row in the matrix as a stepped line (type='S')
#'             
#' @param traces Number of traces to draw if type="traces"
#'             
#' @param new Create a new set of axes (the skyline is fitted to the plotting devie)
#' 
#' @param add Add to the current plot (do not create a new plotting device)
#' 
#' @param ... Parameters passed to plotting function
#' 
#' @export
plotSkyline <- function(times, skyline_mat, type="smooth", traces=1000, col=pal.dark(cblack), fill=pal.dark(cgray, 0.25), 
                        new=TRUE, add=FALSE, xlims=NULL, ylims=NULL, ...) {
  
  
  # Check dimensions
  if (length(times) != ncol(skyline_mat)) stop("Dimension mismatch between times and skyline_mat!")
  
  # Do not open a new device (add to the current plot)
  if (add == TRUE) par(new=TRUE)
  
  # Create a new set of axes that fits this skyline
  if (new == TRUE) {
      if (is.null(ylims)) ylims <- paddedrange(skyline_mat)
      if (is.null(xlims)) xlims <- range(times)
      
      plot(1, type='n', ylim=ylims, xlim=xlims, ...)
  }
  
  #######################
  # Plot actual skyline #
  #######################
  
  # Plot traces
  if (type == "lines" || type == "steplines") {
      if (traces > 0 && traces < nrow(skyline_mat)) 
          ind <- sample(nrow(skyline_mat), size=traces, replace=FALSE)
      else 
          ind <- 1:nrow(skyline_mat)
    
      for (i in 1:length(ind)) {
        lines(times, skyline_mat[ind[i],], col=col, type=if(type=="lines") 'l' else 'S')
      }   
  } else
  # Plot HPDs smooth
  if (type == "smooth" && nrow(skyline_mat) == 3) {
    
      polygon(c(times, rev(times)), c(skyline_mat[1,], rev(skyline_mat[3,])), col=fill, border=NA)
      lines(times, skyline_mat[2,], col=col)  
      
  } else
  # Plot HPDs stepped
  if (type == "step" && nrow(skyline_mat) == 3) {
    
      for (i in 2:ncol(skyline_mat)) {
        rect(times[i-1], skyline_mat[1,i-1], times[i], skyline_mat[3,i-1], col=fill, border=NA)
      }
      lines(times, skyline_mat[2,], col=col, type='S')
      
  } else
    stop("Invalid type parameter for Skyline plot!")
  
  
  
  if (add == TRUE) par(new=FALSE)
}





#' Function to plot a prettier skyline with a lot of options.
#'
#' When plotting the y-axis on the right need to ensure that the plot margins are big enough!
#'
#' @param side Side to draw the y-axis
#' @param xline Line to draw x-axis label on
#' @param yline Line to draw y-axis label on
#' 
#' @export
plotSkylinePretty <- function(times, skyline_mat, type="smooth", traces=1000, col=pal.dark(cblack), fill=pal.dark(cgray, 0.25), col.axis=NA,  
                              xaxis=TRUE, yaxis=TRUE, xlab="", ylab="", xline=1, yline=1, xticks=NULL, yticks=NULL, xlims=NULL, ylims=NULL, axispadding=0, side=2, ...) {
  
  # Plot the basic skyline first
  plotSkyline(times, skyline_mat, type=type, traces=traces, col=col, fill=fill, axes=FALSE, xlab=NA, ylab=NA, xlims=xlims, ylims=ylims,  ...)

  
  # General settings
  if (is.null(xticks)) xticks <- axTicks(1)
  if (is.null(yticks)) yticks <- axTicks(2)
  if (is.null(ylims))  ylims  <- paddedrange(skyline_mat)
  if (is.null(xlims))  xlims  <- range(times)

  
  # y-axis
  if (yaxis == TRUE) {
      if (side == 2) 
        xpos <- min(xlims, xticks)-(axispadding*diff(xlims))
      else
        xpos <- max(xlims, xticks)+(axispadding*diff(xlims))
    
      if (is.na(col.axis)) col.axis <- pal.dark(cblack)
      
      axis(side, at=yticks, pos=xpos, lwd=0, lwd.ticks=1, las=1, col=col.axis, col.axis=col.axis)
      lines(rep(xpos,2), range(ylims, yticks), col=col.axis)
      mtext(ylab, side=side, line=yline, col=col.axis)
  }
  
  
  # x-axis
  if (xaxis == TRUE) {
    ypos <- min(ylims, yticks)-axispadding*diff(ylims)
    axis(1, at=xticks, pos=ypos, lwd=0, lwd.ticks=1, las=1)
    lines(range(xlims, xticks), rep(ypos,2))
    mtext(xlab, side=1, line=xline)
  }  
}


