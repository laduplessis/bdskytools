###############################################################################
# Initialize palettes for plotting

require(RColorBrewer)

cblue   <- 1
cgreen  <- 2
cred    <- 3
corange <- 4
cpurple <- 5
cgray   <- cgrey <- 7
cblack  <- 8
cwhite  <- 9

pal.light <- c(brewer.pal(12,"Paired")[seq(1,12,2)], "#dddddd", "#000000", "#ffffff")
pal.dark  <- c(brewer.pal(12,"Paired")[seq(2,12,2)], "#777777", "#000000", "#ffffff")

pal.trans <- function(col, alpha=0.25) {
  return (paste0(pal.dark[col],format(as.hexmode(round(alpha*255)), width=2)))
}


