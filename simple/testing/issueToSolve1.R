## in the following dataset, even the actual tree has a posterior of zero

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

load("Robjects/data1.RData")


## run outbreaker with only a few iterations
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w.dens+.001, init.tree=dat$ances, n.iter=1)

res$chain$post

plot(res$chains$step, res$chain$post, type="l", main="posterior")

plot(res$chains$step[res$chains$step>2e5], res$chain$post[res$chains$step>2e5], type="l", main="posterior")
