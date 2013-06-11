## in the following dataset, even the actual tree has a posterior of zero

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

## used the make the data
w <- c(0,.2,.3,1,.8,.2,.1,.05)
full <- simOutbreak(R0=2, infec.curve=w)
dat <- full[1:20]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
#save(w, full, dat, collecDates, file="Robjects/data2big.RData")


##load("Robjects/data2.RData")


## run outbreaker, init with seqTrack
system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree=dat$ances, n.iter=5e5))

res$chain$post

plot(res$chains$step, res$chain$post, type="l", main="posterior")

plot(res$chains$step[res$chains$step>1e5], res$chain$post[res$chains$step>1e5], type="l", main="posterior")



## check results ##
plot(density(res$chains$mu1[res$chains$step>2e5]), main="mu1")
plot(density(res$chains$mu2[res$chains$step>2e5]), main="mu2")
x <- get.TTree.simple(res)
plot(x)

## proportion of correct inferrences
mean(x$ances==dat$ances,na.rm=TRUE)
