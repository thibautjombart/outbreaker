## in the following dataset, even the actual tree has a posterior of zero

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

load("Robjects/data3.RData")


## init with seqTrack and random collection dates, actually try to reconstruct outbreak
res <- outbreaker(dna=dat$dna, dates=dat$dates+3, w.dens=w.dens, init.tree=dat$ances, n.iter=3e4)


res$chain$post

## should not go up as initial result is the real tree
plot(res$chains$step, res$chain$post, type="l", main="posterior")

plot(res$chains$step[res$chains$step>2e5], res$chain$post[res$chains$step>2e5], type="l", main="posterior")



## check results ##
plot(density(res$chains$mu1[res$chains$step>2e5]), main="mu1")
plot(density(res$chains$mu2[res$chains$step>2e5]), main="mu2")
x <- get.TTree.simple(res)
plot(x)
