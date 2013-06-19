
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)
## note: this works only if outbreak has at least 30 case
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100)[1:15]
collecDates <- dat$dates + sample(0:3, size=length(dat$dates), replace=TRUE, prob=w)
#plot(dat)
D <- prop.table(matrix(1:length(dat$dates)^2, ncol=length(dat$dates), nrow=length(dat$dates)),1)


## test parallel
#res <-  outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5)


## test spatial matrix
#res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5)

res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=5e4, init.spa1=10, init.spa2=69.69, spa.model=1)

dat$ances
get.tTree(res)$ances
