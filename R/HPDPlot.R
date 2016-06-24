###############################################################################
# Functions for plotting sets of HPDs next to each other


#' Plot a rounded bar (jelly bean)
#'
#' @param x  x-coordinate of the bar
#' @param y1 Upper coordinate of the bar
#' @param y2 Lower coordinate of the bar
#' @param width Width of the bar
#' @param nv Number of vertices to use for drawing the circle
#' @param plotAlways Plot even if the proportions will make it look bad
#' @param flip Plot horizontally (x is interpreted as y, and y1, y2 as x1, x2)
#' @param side Plot only half of the bar. (If side = 2 plot the left half, if side = 4 plot the right half)
#' 
#' @export
plotRoundedBars <- function(x, y1, y2, width, nv=100, plotAlways=FALSE, flip=FALSE, side=NA, ...) {
    
  pin     <- par("pin")
  usrdiff <- diff(par("usr"))
  scale   <- pin[1]/pin[2] * usrdiff[3]/usrdiff[1]
  if (flip == TRUE) 
    scale <- 1/scale

  radius <- width/2
  
  angles <- seq(0, pi, length.out=nv)    
  if (!is.na(side) && !flip) {
    if (side == 4) {
      angles <- seq(0, pi/2, length.out=nv/2)    
    } else 
      if (side == 2) {
        angles <- seq(pi/2, pi, length.out=nv/2)        
    } 
  }
      
  for (i in 1:length(x)) {
    xv    <- cos(angles)*radius + x[i]
    xv    <- c(xv, rev(xv))    
    ytemp <- sin(angles)*radius*scale
    
    if (!is.na(side) && flip) {
      if (side == 2) { 
        yv <- c(rep(y2[i]-0.5*abs(y1[i]-y2[i]), length(ytemp)), y1[i]+(radius*scale)-rev(ytemp))          
      } else
        if (side == 4) {        
          yv <- c(y2[i]-(radius*scale)+ytemp, rep(y1[i]+0.5*abs(y1[i]-y2[i]),length(ytemp)))                      
        }              
    } else {
      yv <- c(y2[i]-(radius*scale)+ytemp, y1[i]+(radius*scale)-rev(ytemp))          
    }
        
    if (plotAlways || width*scale < abs(y2[i]-y1[i])) {                                  
      if (flip == TRUE) {
        polygon(yv, xv, ...)
      } else 
        polygon(xv, yv, ...)
    }
  }
}


getBeanHPDs <- function(data) {
  
  hpd1 <- getMatrixHPD(data, alpha=0.5)
  hpd2 <- getMatrixHPD(data, alpha=0.05)
  
  y <- rbind(hpd2[1,], hpd1, hpd2[3,])
  rownames(y) <- c(2.5, 25, 50, 75, 97.5)
  
  return(y)
}


#' Plot jelly beans
#'
#'
#' @param x locations to plot jellybeans
#' @param y matrix with 5 rows and n columns
#'                    
#' @export
plotJellyBeans <- function(x, y, maxwidth=0.5, plotmedian=TRUE, plotinterquartile=TRUE, col=pal.dark(cblue), border=NA, side=NA) {
  
  y <- as.matrix(y)
  
  plotRoundedBars(x, y[1,], y[5,], width=0.5*maxwidth, col=col, border=border, side=side)
  
  if (plotinterquartile == TRUE) {
    plotRoundedBars(x, y[2,], y[4,], width=maxwidth, col=col, border=border, side=side)
  }
  
  if (plotmedian == TRUE) { 
    midwidth <- min(y[3,]-y[1,], y[5,]-y[3,], diff(par("usr"))[3]*0.025) 
    plotRoundedBars(y[3,], x-0.75*maxwidth, x+0.75*maxwidth, width=midwidth, col=pal.dark(cwhite), border=border, flip=TRUE, plotAlways=TRUE, side=side)  
    plotRoundedBars(y[3,], x-0.375*maxwidth, x+0.375*maxwidth, width=midwidth/3, col=col, border=border, flip=TRUE, plotAlways=TRUE, side=side)    
  }
}


#' Get HPDs and plot jelly beans
#'
#'
#' @param data       List with each element being a vector that represents a posterior sample
#'                   Or a matrix, with each column being a  posterior sample
#'                    
#' @export
plotJellyBeanHPDs <- function(data, maxwidth=0.5, hpdonly=FALSE, plotmedian=TRUE, plotinterquartile=TRUE, col=pal.dark(cblue), border=NA, side=NA) {
  
  x <- 1:ncol(data)
  y <- getBeanHPDs(data)
  
  if (hpdonly == FALSE) {
      datarange <- apply(data, 2, range)
      segments(1:ncol(data), datarange[1,], 1:ncol(data), datarange[2,], col=border)
  }
  
  plotJellyBeans(x,y, maxwidth=maxwidth, plotmedian=plotmedian, plotinterquartile=plotinterquartile, col=col, border=border, side=side)
}




#' Plot beans
#' 
#' @param data       List with each element being a vector that represents a posterior sample
#'                   Or a matrix, with each column being a  posterior sample
#'                   
#' @export
plotBeanPlot <- function(data, hpdlines=c(0.05, 0.5), linewidths=c(0.5,0.75), maxwidth=0.5, hpdonly=TRUE, plotmedian=TRUE, 
                         col=pal.dark(cgray, 0.25), border=pal.dark(black), side, ...) {
  
    
  beanside <- 'no'
  offset <- c(-0.5,0.5)
  if (!is.na(side)) {
    if (side == 2) {
      beanside <- 'first'
      offset <- c(-0.5,0)
    } else {
      beanside <- 'second' 
      offset <- c(0.5,0)
    }
  }
  
  beanplot::beanplot(data, what=c(0,1,0,0), col=col, border=border, add=TRUE, frame.plot=FALSE, axes=FALSE, width=maxwidth, side=beanside, ...)

  # Plot median
  if (plotmedian == TRUE) {
    width <- 1.25*maxwidth      
    for (i in 1:ncol(data)) {
      med <- median(data[[i]])
      lines(i+offset*width, rep(med,2), col=border)
      lines(i+offset*width, rep(med,2), col=border)            
    }
  }
  
  # Plot lines
  for (i in 1:length(hpdlines)) {
      width <- maxwidth*linewidths[i]
      hpd   <- getMatrixHPD(data, alpha=hpdlines[i])
      segments(1:ncol(data)+offset[1]*width, hpd[1,], 1:ncol(data)+offset[2]*width, hpd[1,], col=border)
      segments(1:ncol(data)+offset[1]*width, hpd[3,], 1:ncol(data)+offset[2]*width, hpd[3,], col=border)    
  }
  
}



#' Plot distributions using beans
#'
#' @param data       List with each element being a vector that represents a posterior sample
#'                   Or a matrix, with each column being a  posterior sample
#' @param type       One of 'simple', 'beans' or 'jellybeans'       
#' @param hpdlines   Alpha values for HPDs to draw horizontal lines for (does not apply to type='jellybeans')       
#' @param hpdonly    Only plot the biggest HPD in HPD lines
#' @param showlimits Add the upper/lower limits to the margin for HPDs out of bounds
#' @param maxwidth   The maximum width of the HPD (does not apply to type='simple')
#' 
#'
#' @export
plotHPDs <- function(data, type='beans', hpdlines=c(0.05, 0.5), linewidths=c(0.75,1.0), maxwidth=0.5, hpdonly=TRUE, showlimits=FALSE, plotmedian=TRUE, plothalf=FALSE,
                     col=pal.dark(cblack), fill=pal.dark(cgray, 0.25), col.axis=NA, lwd=1, add=FALSE, new=TRUE,
                     xaxis=TRUE, yaxis=TRUE, xlab="", ylab="", xline=1, yline=1, xticks=NULL, yticks=NULL, xticklabels=NULL, 
                     ygrid=FALSE, xgrid=FALSE, gridcol=pal.light(cgray), 
                     ylims=NULL, cex.label=1.4, cex.axis=1, axispadding=0, side=2, ...) {
  
  
    ##################
    # Initialization #
    ##################
  
    xlims <- c(0, ncol(data)+1)
        
    # Do not open a new device (add to the current plot)
    if (add == TRUE) par(new=TRUE)
    
    # Create a new set of axes    
    if (new == TRUE) {      
        if (is.null(ylims)) ylims <- paddedrange(data)
      
        plot(1,type="n",xlim=xlims,ylim=ylims, axes=FALSE, ylab=NA, xlab=NA)                     
    } else {
        if (is.null(ylims)) ylims <- range(axTicks(2))      
    }

    # Get and set clipping area
    usr <- par("usr")
    clip(xlims[1], xlims[2], ylims[1], ylims[2])        

    
    #########################
    # Plot the actual beans #
    #########################
    
    if (type == "jellybeans")
        plotJellyBeanHPDs(data, hpdonly=hpdonly, plotmedian=plotmedian, plotinterquartile=0.5 %in% hpdlines, col=fill, maxwidth=maxwidth, border=col, side=if (plothalf) side else NA)
    else 
    if (type == "beans")
        plotBeanPlot(data, hpdlines=hpdlines, linewidths=linewidths, maxwidth=maxwidth, hpdonly=hpdonly, plotMedian=plotMedian, 
                     col=fill, border=col, side=if (plothalf) side else NA, ...)    
    else 
        stop("Invalid type parameter for HPD plot!")
        
    #######################     
    # Plot axes and grids #
    #######################
    
    # General settings    
    if (is.null(xticks)) xticks <- axTicks(1)
    if (is.null(yticks)) yticks <- axTicks(2)
    ypos  <- min(ylims, yticks)-axispadding*diff(ylims)
    xpos2 <- min(xlims, xticks)-(axispadding*diff(xlims))
    xpos4 <- max(xlims, xticks)+(axispadding*diff(xlims))
        
    # x-grid
    if (xgrid == TRUE) {
      for (tick in yticks) lines(c(xpos2,xpos4), rep(tick,2), col=gridcol, lty=3, lwd=0.5)
    }
    
    # y-grid
    if (ygrid == TRUE) {
      for (x in 1:ncol(data)) lines(rep(x,2), c(ypos,max(yticks)), col=gridcol, lty=3, lwd=0.5)      
    }
    
    # x-axis
    if (xaxis == TRUE) { 
      if (is.null(xticklabels)) 
        xticklabels <- TRUE
      
      axis(1, at=xticks, labels=xticklabels, pos=ypos, lwd=0, lwd.ticks=lwd, las=1, cex.axis=cex.axis)
      lines(range(xlims, xticks), rep(ypos,2), lwd=lwd)
    }
    mtext(xlab, side=1, line=xline, cex=cex.label)
        
    # y-axis
    if (yaxis == TRUE) {
      if (side == 2) 
        xpos <- xpos2
      else
        xpos <- xpos4
      
      if (is.na(col.axis)) col.axis <- pal.dark(cblack)
      
      axis(side, at=yticks, pos=xpos, lwd=0, lwd.ticks=lwd, las=1, col=col.axis, col.axis=col.axis, cex.axis=cex.axis)
      lines(rep(xpos,2), range(ylims, yticks), col=col.axis, lwd=lwd)
    }
    mtext(ylab, side=side, line=yline, col=col.axis, cex=cex.label)
    
    # Plot upper and lower limits that fall outside axis
    if (showlimits == TRUE) {

        at <- 1:ncol(data)  
        if (plothalf == TRUE) 
          at <- at + (side-3)/2
              
        upper  <- apply(data, 2, max)
        labels <- sapply(upper, function(x) sprintf("%.2f",x))
        labels[which(upper <= ylims[2])] = ""
        axis(3, at=at, labels=labels, las=2, lwd=0, lwd.ticks=0, pos=ylims[2], col.axis=pal.dark(cred), cex.axis=0.5*cex.axis)
        
        lower  <- apply(data, 2, min)
        labels <- sapply(lower, function(x) sprintf("%.2f",x))
        labels[which(lower >= ylims[1])] = ""
        axis(1, at=at, labels=labels, las=2, lwd=0, lwd.ticks=0, pos=ylims[1], col.axis=pal.dark(cred), cex.axis=0.5*cex.axis)
        
        #labels <- c()
        #for (m in maxes) {
        #    if (m > ylims[2]) 
        #        labels <- c(labels, sprintf("%.2f",m)) 
        #    else
        #        labels <- c(labels, "")        
        #}                
    }    
    
    ##################
    # Reset settings #
    ##################
    
    if (add == TRUE) par(new=FALSE)
    do.call("clip", as.list(usr))
        
}
