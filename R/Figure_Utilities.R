###############################################################################
# Initialize palettes for plotting

#requireNamespace(RColorBrewer)

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
cgray   <- cgrey <- 7

#' @export
cblack  <- 8

#' @export
cwhite  <- 9

palette.dark  <- c(RColorBrewer::brewer.pal(12,"Paired")[seq(2,12,2)], "#777777", "#000000", "#ffffff")
palette.light <- c(RColorBrewer::brewer.pal(12,"Paired")[seq(1,12,2)], "#dddddd", "#000000", "#ffffff")

#' @export
pal.light  <- function(col, alpha=1) {
  return (paste0(palette.dark[col],format(as.hexmode(round(alpha*255)), width=2)))
}

#' @export
pal.dark  <- function(col, alpha=1) {
  return (paste0(palette.dark[col],format(as.hexmode(round(alpha*255)), width=2)))
}



