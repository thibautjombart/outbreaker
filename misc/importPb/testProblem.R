
library(outbreaker)
#load("bugSmall.RData")

load("bug.RData")
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, idx.dna=which.seq, n.iter=5e4)
