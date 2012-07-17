
## load packages ##
library(outbreaker)
library(adegenet)

## get simulated outbreak, >=15 cases
dat <- list(n=0)

## curve of infectiousness from t=0 to...
w.dens <- c(0,1,2,1)

while(dat$n<15){
    dat <- simOutbreak(R0=2, infec.curve=w.dens, n.hosts=30,tree=FALSE)
}

## check outbreak
dat

## plot
## code vertices: red=old, white=recent
## code edges: blue = not many mutations, red=largest nb of mutations
plot(dat, main="simulated outbreak")

## make collection dates
collecDates <- dat$dates + sample(0:3, dat$n, prob=w.dens, replace=TRUE)


## try to reconstruct outbreak using outbreaker ##
## ! THIS TAKES TIME ! ##
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w.dens+.001, init.tree="seqTrack", n.iter=3e5)


## check that we got off zero likelihood
res$chain$post
plot(res$chains$step, res$chain$post, type="l", main="posterior")
plot(res$chains$step[res$chains$step>2e5], res$chain$post[res$chains$step>2e5], type="l", main="posterior")


## get summarized results
x <- get.TTree.simple(res, burnin=2e5) # large burnin


## accuracy of ancestries
mean(x$ances==dat$ances, na.rm=TRUE) # proportion of correct ancestries
weighted.mean(x$ances==dat$ances, na.rm=TRUE,w=x$p.ances) # same, weighted by posterior proba (should be same or better)
boxplot(x$p.ances ~ x$ances==dat$ances, main="posterior proba for ancestries", xlab="accurate ancestry", ylab="posterior probability")

## nb of generations
x$n.gen # should be only 1s

## mutation rates
par(mfrow=c(2,1))
plot(density(res$chains$mu1[res$chains$step>2e5]), main="mu1 post")
abline(v=1e-4)
mtext(side=3, paste("actual value:",1e-4))

plot(density(res$chains$mu2[res$chains$step>2e5]), main="mu2 post")
abline(v=0.5e-4)
mtext(side=3, paste("actual value:",0.5e-4))


## infection times
Tinf <- res$chains[res$chains$step>2e5,grep("Tinf",names(res$chains))]
library(UsingR)
violinplot(Tinf)
points(dat$dates, col="red")


## plot reconstructed outbreak
plot(x)

## graph-like plot
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(dat,main="simulated data")
plot(g <- as.igraph(x), main="reconstructed outbreak")

