library(bdskytools)

fname <- "~/Documents/Projects/Ebola/Data/HCV_OUPrior/HCV_oup_40_127.log"

lf <- readLogfile(fname)

R0_sky <- getSkylineSubset(lf, "R0")
delta_sky <- getSkylineSubset(lf, "becomeUninfectiousRate")

R0_hpd <- getMatrixHPD(R0_sky)


timegrid      <- 1:500
R0_gridded    <- gridSkyline(R0_sky, lf$orig_root, timegrid)
delta_gridded <- gridSkyline(delta_sky, lf$orig_root, timegrid)

R0_gridded_hpd    <- getMatrixHPD(R0_gridded)
delta_gridded_hpd <- getMatrixHPD(delta_gridded)


plotSkylinePretty(timegrid, R0_gridded_hpd, axispadding=0.0, ylims=c(0,2), col=pal.dark(corange), fill=pal.dark(corange,0.25), col.axis=pal.dark(corange), xlab="Time", ylab="R0", yline=2, xline=2)
plotSkylinePretty(timegrid, delta_gridded_hpd, axispadding=0.0, col=pal.dark(cpurple), fill=pal.dark(cpurple,0.25), col.axis=pal.dark(cpurple), xaxis=FALSE, ylab="Delta", side=4, yline=2, add=TRUE)



plotSkyline(timegrid, R0_gridded_hpd, type='smooth')
plotSkyline(timegrid, R0_gridded_hpd, type='step')
plotSkyline(1:40, R0_hpd, type='steplines')

x <- c(seq(0,10,2), seq(20, 40, 3), 40:66)
plotSkyline(c(seq(0,10,2), seq(20, 40, 3), 40:66), R0_hpd, type='step')
for (y in x) abline(v = y)

plotSkyline(timegrid, R0_gridded_hpd, type='smooth', col=pal.dark(cblue), fill=pal.dark(cblue,0.25))
plotSkyline(timegrid, R0_gridded, type='steplines', col=pal.dark(cgreen,0.1), traces=100, add=TRUE, new=FALSE)

plotSkyline(1:40, revMatrix(R0_sky[1:3,]), type='steplines', col=pal.dark(cgreen,1), traces=3)
plotSkyline(timegrid, R0_gridded[1:3,], type='steplines', col=pal.dark(cblue,1), traces=3, add=TRUE, new=TRUE)



plotSkylinePretty(timegrid, R0_gridded_hpd, axispadding=0.0, ylims=c(0,2), col=cblue, col.axis=pal.dark(cblue), xlab="time", ylab="R0", yline=2, xline=2)

