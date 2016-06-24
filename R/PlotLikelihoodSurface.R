#library(plot3D)



#' Remove all rows containing non-finite values
removeInfinite <- function(data) {
  idxs <- c()
  for (i in 1:ncol(data))
    idxs <- union(idxs, which(is.finite(data[,i]) == FALSE))
  print(idxs)
  print(data$popSize[idxs])
  
  if (length(idxs) > 0) {
    print(paste("Removed", length(idxs), "rows"))
    return(data[-idxs,])
  } else
    return(data)
}


#' Return datastructure with the meshed likelihood surface as well as the maximum
#' 
#' Does not assume par1 and par2 are ordered (so the code is not very efficient)
#' 
#' @export
getSurface <- function(xpar, ypar, lk, crop=TRUE) {
  
  x  <- sort(unique(xpar))
  y  <- sort(unique(ypar))
  xy <- plot3D::mesh(x, y)
  z  <- matrix(nrow=length(x), ncol=length(y))
  
  maxlk <- -Inf
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      z[i,j] <- lk[which(xpar == x[i] & ypar == y[j])]
      if (is.finite(z[i,j]) && z[i,j] > maxlk) {
          maxlk <- z[i,j]
          maxx  <- i
          maxy  <- j
      }
    }
  }

  max <- c(x[maxx], y[maxy], z[maxx,maxy])
  
  if (crop == TRUE) {
      i <- length(x)
      while (sum(is.finite(z[i,])) == 0) {        
        i <- i - 1
      }
      
      j <- length(y)
      while (sum(is.finite(z[,j])) == 0) {
        j <- j - 1
      }
    
      #print(paste("Cropped from ",length(x),"x",length(y)," to ",i,"by",j))
      return(list(x=x[1:i], y=y[1:j], z=z[1:i,1:j], xy=xy, max=max))
  } else
      return(list(x=x, y=y, z=z, xy=xy, max=max))

}




#' Plot cross section of the likelihood surface of par1, for given values of par2
#' Logfile is a BEAST2 log file
#' 
#' @export
plotCrossSection <- function(par1, par2, lk, lk.upper=NULL, lk.lower=NULL, par1.truth=NA, par2.truth=NA,
                             col=pal.dark(cblue), fill=pal.dark(cblue,0.5),  
                             xlims=NULL, ylims=NULL, par1.label=NA, par2.label=NA, lk.label=NA, sections=-1, ...) {
  
  plotsymbols   <- c(15,0, 16 , 1, 17,2, 18, 6, 3, 4, 8, 9, 10, 13)
  
  #if (!is.expression(par2.label && is.na(par2.label)) 
  #  par2.label = ''
  
  if (is.null(xlims))
    xlims <- range(pretty(par1))
  
  if (is.null(ylims))
    ylims <- range(pretty(lk))
    
  if (length(sections) > 1) {
      par2.unique  <- sections
      plotsections <- seq(1,length(sections))
  } else {
      par2.unique <- unique(par2)
  
      if (sections > 0) 
          plotsections <- seq(1,length(par2.unique), length.out=min(length(par2.unique), sections))
      else
          plotsections <- seq(1,length(par2.unique), length.out=length(par2.unique))
  }
  
  
  plot(1,type='n', xlim=xlims, ylim=ylims, bty='n', xlab=par1.label, ylab=lk.label, ...)
  grid(col='lightgrey')
  
  legend   = c()
  legend.y = c()

  for (i in plotsections) {
    idxs <- which(par2 == par2.unique[i])
  
    if (par2.unique[i] == par2.truth) {
        plotcol  <- pal.dark(cred)
        fillcol  <- pal.dark(cred,0.5)  
        #legend.y <- lk[idxs[which(par1[idxs] == max(par1[idxs]))]]
    } else {
        plotcol <- col
        fillcol <- fill
    }
    
    if (!is.null(lk.upper) && !is.null(lk.lower)) {          
        arrows(par1[idxs], lk.lower[idxs], par1[idxs], lk.upper[idxs], length=0.01, angle=90, code=3, col=fillcol)
    }
    
    points(par1[idxs], lk[idxs], col=plotcol, pch=plotsymbols[i %% length(plotsymbols)])
    lines(par1[idxs], lk[idxs], col=fillcol, lwd=1)
    
    
    
    #legend.y <- c(legend.y, lk[idxs[which(par1[idxs] == max(par1[idxs]))]])
    #legend   <- c(legend, paste0(par2.label," = ", par2.unique[i]))
    #legend <- c(legend, bquote(.(par2.label[[1]])~" = 5"))
  }
  mtext(par2.label,side=4,line=1)
  #axis(4, at=legend.y, labels=par2.label, tick=FALSE, lwd=0, las=1, line=-1)
  #if (sum(is.finite(legend.y)) > 0)
  #    axis(4, at=legend.y, labels=eval(legend), tick=FALSE, lwd=0, las=1, line=-1)
    
  
  abline(v=par1.truth, col=pal.dark(cred), lty=2)
  
  
}

