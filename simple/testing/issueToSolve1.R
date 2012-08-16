## in the following dataset, even the actual tree has a posterior of zero

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

load("Robjects/data1.RData")


## init with actual tree and non-random, most likely collection dates
res <- outbreaker(dna=dat$dna, dates=dat$dates+2, w.dens=w.dens+.001, init.tree=dat$ances, n.iter=5)


## init with actual tree and random collection dates
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w.dens+.001, init.tree=dat$ances, n.iter=5)


## init with seqTrack and random collection dates, actually try to reconstruct outbreak
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w.dens+.001, init.tree="seqTrack", n.iter=2e6)


res$chain$post

plot(res$chains$step, res$chain$post, type="l", main="posterior")

plot(res$chains$step[res$chains$step>2e5], res$chain$post[res$chains$step>2e5], type="l", main="posterior")



## check results ##
plot(density(res$chains$mu1[res$chains$step>2e5]), main="mu1")
plot(density(res$chains$mu2[res$chains$step>2e5]), main="mu2")
x <- get.TTree.simple(res)
plot(x)
