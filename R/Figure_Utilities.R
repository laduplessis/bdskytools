###############################################################################
# Initialize colour constants

#' @export
cblue   <- 1

#' @export
cgreen  <- 2

#' @export
cred    <- 3

#' @export
corange <- 4

#' @export
cpurple <- 5

#' @export
cgray   <- 7

#' @export
cgrey <- 7

#' @export
cblack  <- 8

#' @export
cwhite  <- 9


###############################################################################
# Initialize palettes for plotting

palette.dark  <- c(RColorBrewer::brewer.pal(12,"Paired")[seq(2,12,2)], "#777777", "#000000", "#ffffff")
palette.light <- c(RColorBrewer::brewer.pal(12,"Paired")[seq(1,12,2)], "#dddddd", "#000000", "#ffffff")

#' Palette of diverging light colours (Paired palette from RColorBrewer)
#' 
#' @export
pal.light  <- function(col, alpha=1) {
  return (paste0(palette.light[col],format(as.hexmode(round(alpha*255)), width=2)))
}

#' Palette of diverging dark colours (Paired palette from RColorBrewer)
#' 
#' @export
pal.dark  <- function(col, alpha=1) {
  return (paste0(palette.dark[col],format(as.hexmode(round(alpha*255)), width=2)))
}


###############################################################################
# Plot PDF figures

#' Need to follow with dev.off()
#' 
#' @export
NewFig <- function(figname="fig.pdf", width=3.42, aspectratio=3/2, pointsize=6) {
  
  pdf(file=figname, width=width, height=width/aspectratio, pointsize=pointsize, useDingbats=FALSE)
  
  # Make figure
  # dev.off()
}


#' Label plot in top left corner
#' 
#' @export
labelPlot <- function(label, cex=1) {
  
  mfrow <- par("mfrow")
  scale <- mfrow[2]
  
  mar   <- par("mar")
  fin   <- par("fin")
  mai   <- par("mai")
  
  plotheight <- (fin[2]/scale)
  #plotheight <- par("pin")[2]
  textheight <- strheight(label, cex=1, units="inches")
  
  lines <- (plotheight/textheight)-2
  
  #print(paste(plotheight, textheight, lines, par("lheight")))
  
  mtext(label, 2, adj=0, las=1, padj=-(lines/2), line=mar[2], cex=cex)
  
  #for (i in -20:20)
  #  mtext(i, 2, adj=0, las=1, padj=i, line=2, cex=1.5)
}