## TRY TO SEE HOW THE METHOD BEHAVES ##


## CASE 1: NO GENETIC INFORMATION, ALL INFECTIONS AFTER 1 DAY ##

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

## generate simple data
w <- c(0,1)
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=0)
dat <- full[1:20]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)

## init with seqTrack and random collection dates, actually try to reconstruct outbreak
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree="none", n.iter=5e5)

res$chain$post

## should not go up as initial result is the real tree
plot(res$chains$step, res$chain$post, type="l", main="posterior")

plot(res$chains$step[res$chains$step>2e5], res$chain$post[res$chains$step>2e5], type="l", main="posterior")



## check results ##
par(mfrow=c(2,1))
plot(density(res$chains$mu1[res$chains$step>2e5]), main="mu1")
plot(density(res$chains$mu2[res$chains$step>2e5]), main="mu2")

x <- get.TTree.simple(res)

mean(x$ances==dat$ances,na.rm=TRUE)

## support should be 1/number of cases at time step before
par(mfrow=c(1,1))
plot(dat)
x11()
plot(x, annot="prob") # all right thus far
