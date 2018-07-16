library(bdskytools) 

###############################################################################
# Plot theoretical quantiles of a bunch of OU-processes, varying nu and sigma

# Time intervals
dt <- 0.005
tn <- 1
t  <- seq(from=0, to=tn, by=dt)

# Plot processes
n    <- 4
ylim <- c(0,5)
par(mfrow=c(n,n))
for (sigma in seq(0.25,1,length.out=n)) {
    for (nu in seq(0.2,5,length.out=n)) {
        title=substitute(x[0]~"=3,"~mu~"=1, "~sigma~"="~s~", "~nu~"="~n,list(s=sigma,n=nu))
        plotOUProcessHPD(x0=3, t, mu=1, sigma, nu, ylim=ylim, main=title, xlab="t",ylab="x")      
    }
}

