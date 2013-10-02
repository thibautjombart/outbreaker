library(outbreaker)

load("837247b918.RData")

## res <- outbreaker(dat$dna[1:30, ], collecDates[1:30], mut.model=2, w.dens=res$w,
##                  n.iter=5e4, seed=2, find.import=FALSE)

res <- outbreaker(dat$dna, collecDates, mut.model=2, w.dens=res$w)

res <- outbreaker(dat$dna, collecDates, mut.model=2, w.dens=res$w,
                  n.iter=5e4, seed=1, find.import=FALSE)
