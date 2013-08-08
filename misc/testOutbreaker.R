
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)
set.seed(2)

## NON-SPATIAL SIMULATION ##

## this may generate an error if outbreak doesn't take off
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100, spatial=FALSE, mu.transi=2e-4)[1:30]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)


## test parallel
#res <-  outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5)

## run outbreaker
res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4, init.spa1=1, init.spa2=10, spa.model=0)
res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4, init.spa1=1, init.spa2=10, spa.model=0, find.import=FALSE)

dat$ances
get.tTree(res)$ances
mean(dat$ances == get.tTree(res)$ances, na.rm=TRUE)



## SPATIAL SIMULATION ##
## test spatial model 1 (exponential)
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)

## this may generate an error if outbreak doesn't take off
dat <- simOutbreak(R0 = 1.5, infec.curve = w, n.hosts = 100, spatial=TRUE)[1:15]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)
D <- as.matrix(dist(dat$xy))

res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=5e4, spa.model=1, find.import=FALSE)

dat$ances
get.tTree(res)$ances
mean(dat$ances == get.tTree(res)$ances, na.rm=TRUE)

dist.inf <- na.omit(sapply(1:dat$n, function(i) D[dat$id[i], dat$ances[i]]))
mean(dist.inf)

par(mfrow=c(2,2))
plotChains(res)
plotChains(res, what="spa1")
plotChains(res, what="spa1", type="dens")
abline(v=mean(dist.inf))
hist(dist.inf , col="grey",
     main="Distance of infectious contacts", xlab="Distance infector-infected")




## bigger example ##
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)

## this may generate an error if outbreak doesn't take off
dat <- simOutbreak(R0 = 1.5, infec.curve = w, n.hosts = 200, spatial=TRUE)[1:50]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)
D <- as.matrix(dist(dat$xy))

res <-  outbreaker.parallel(5, dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5, spa.model=1)


dat$ances
get.tTree(res)$ances
mean(dat$ances == get.tTree(res)$ances, na.rm=TRUE)

dist.inf <- na.omit(sapply(1:dat$n, function(i) D[dat$id[i], dat$ances[i]]))
mean(dist.inf)

par(mfrow=c(2,2))
plotChains(res)
plotChains(res, what="spa1")
plotChains(res, what="spa1", type="dens")
abline(v=mean(dist.inf))
hist(dist.inf , col="grey",
     main="Distance of infectious contacts", xlab="Distance infector-infected")


