## TRY TO SEE HOW THE METHOD BEHAVES ##


## DATA SIMULATION ##

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

## w <- c(0,.1,.2,.5,2,.5,.2,.1)
## full <- simOutbreak(R0=2, infec.curve=w, mu.transi=2e-4, mu.transv=1e-4)
## dat <- full[1:20]
## collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
## save(w, full, dat, collecDates, file="Robjects/data4.RData")

load("Robjects/data4.RData")
plot(dat, main="Data")






## ESTIMATE MUTATION RATES ONLY ##
## run outbreaker
system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree=dat$ances, n.iter=5e5, move.mut = TRUE,
                              move.ances = FALSE, move.kappa = FALSE, move.Tinf = FALSE, move.pi = FALSE, move.phi = FALSE))

## check results
plot(res$chains$step, res$chain$post, type="l", main="posterior")

par(mfrow=c(2,1))
plot(density(res$chains$mu1[res$chains$step>2e5]), main="mu1")
abline(v=2e-4, col="blue")
plot(density(res$chains$mu2[res$chains$step>2e5]), main="mu2")
abline(v=1e-4, col="blue")

x <- get.TTree.simple(res)
mean(x$ances==dat$ances,na.rm=TRUE)






## ESTIMATE ALPHA ONLY ##
## run outbreaker
system.time(res <- outbreaker(dna=dat$dna, dates=dat$dates+4, w.dens=w, init.tree="none", n.iter=5e5, init.mu1=2e-4,
                              init.gamma=0.5, move.mut = FALSE, move.ances = TRUE, move.kappa = FALSE,
                              move.Tinf = FALSE, move.pi = FALSE, move.phi = FALSE))

## check results
plot(res$chains$step, res$chain$post, type="l", main="posterior")

x <- get.TTree.simple(res)
mean(x$ances==dat$ances,na.rm=TRUE)

v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$ances!=dat$ances)
if(length(notOk)>0) v.col[notOk] <- "red"
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction", vertex.color=v.col)







## ESTIMATE KAPPA ONLY (FLAT PRIOR ON PI) ##
## run outbreaker
system.time(res <- outbreaker(dna=dat$dna, dates=dat$dates+4, w.dens=w, init.tree=dat$ances, n.iter=5e5, init.mu1=2e-4,
                              init.gamma=0.5, pi.param1 = 1, pi.param2 = 2, move.mut = FALSE, move.ances = FALSE,
                              move.kappa = TRUE, move.Tinf = FALSE, move.pi = FALSE, move.phi = FALSE))

## check results
plot(res$chains$step, res$chain$post, type="l", main="posterior")

x <- get.TTree.simple(res)
mean(x$n.gen==1,na.rm=TRUE)

v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$n.gen!=1)
if(length(notOk)>0) v.col[notOk] <- "red"
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction", vertex.color=v.col)








## ESTIMATE TINF ONLY ##
## run outbreaker
system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree=dat$ances, n.iter=5e5, init.mu1=2e-4,
                              init.gamma=0.5, move.mut = FALSE, move.ances = FALSE,
                              move.kappa = FALSE, move.Tinf = TRUE, move.pi = FALSE, move.phi = FALSE))

## check results
plot(res$chains$step, res$chain$post, type="l", main="posterior")

toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>1e5, toKeep]
points(dat$dates, col="red", pch="x",cex=2)





