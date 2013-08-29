
library(outbreaker)

##load("bugSmall.RData")

load("bug.RData")

which.seq <- c(1,sample(2:114, 10))
which.seq <- c(1,1)
res <- outbreaker(dna=dat$dna[which.seq], dates=collecDates, w.dens=w, idx.dna=which.seq, n.iter=5e4)
