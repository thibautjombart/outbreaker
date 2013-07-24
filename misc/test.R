
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)


## NON-SPATIAL SIMULATION ##
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100, spatial=FALSE)[1:15]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)


## test parallel
#res <-  outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5)

## run outbreaker
res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4, init.spa1=1, init.spa2=10, spa.model=0)

dat$ances
get.tTree(res)$ances



## test spatial model 1 (exponential)
dat <- simOutbreak(R0 = .5, infec.curve = w, n.hosts = 100, spatial=TRUE)[1:15]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)
D <- as.matrix(dist(dat$xy))

res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=5e4, init.spa1=10, init.spa2=69.69, spa.model=1)
