
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)
## note: this works only if outbreak has at least 30 case
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100)[1:15]
collecDates <- dat$dates + sample(0:3, size=15, replace=TRUE, prob=w)
plot(dat)

res <-  outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates,w.dens=w, n.iter=1e5)

dat$ances
get.tTree(res)$ances
