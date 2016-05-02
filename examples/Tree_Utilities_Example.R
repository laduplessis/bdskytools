#rm(list = ls())
library(TreeSim)

plotTreeTimes <- function(tree, height, times, col) {
  
  plot(1,type='n', xlim=c(0,height), ylim=c(0,tree$Nnode))
  par(new=TRUE)
  plot(tree, root.edge=TRUE)
  
  for (i in 1:length(times)) {
      abline(v = times[i], col=col, lty=2)
  }
  
}

########################
# Present tree example #
########################
set.seed(1)
treepresent   <- sim.bd.taxa(n=10, numbsim=1, lambda=2, mu=0.5, frac=1, complete=FALSE, stochsampling=TRUE)

times <- getTreeIntervals(treepresent[[1]], includeRoot=TRUE)
plotTreeTimes(treepresent[[1]], 5, times[which(times[,2] == 1),1], pal.dark(cblue))
plotTreeTimes(treepresent[[1]], 5, times[which(times[,2] == 0),1], pal.dark(cred))




###########################
# Serial sampling example #
###########################
set.seed(15)
treestt <- sim.bdsky.stt(n=10, lambdasky=c(2,1,2), deathsky=c(1,0.5,1.5), timesky=c(0,1,2), sampprobsky=c(0.5,0.5,0.5), rho=0, timestop=0);

times <- getTreeIntervals(treestt[[1]], includeRoot=TRUE)
plotTreeTimes(treestt[[1]], 5, getBranchingTimes(treestt[[1]]), pal.dark(cblue))
plotTreeTimes(treestt[[1]], 5, getLeafTimes(treestt[[1]]), pal.dark(cred))




######################
# Mixed tree example #
######################
set.seed(1)
treemixed <- sim.bdsky.stt(n=0, lambdasky=c(2,1,2), deathsky=c(1,0.5,1.5), timesky=c(0,1,2), sampprobsky=c(0.5,0.5,0.5), rho=0.5, timestop=2.5);

times <- getTreeIntervals(treemixed[[1]])
plotTreeTimes(treemixed[[1]], 5, times[which(times[,2] == 1),1], pal.dark(cblue))
plotTreeTimes(treemixed[[1]], 5, times[which(times[,2] == ),1], pal.dark(cred))
