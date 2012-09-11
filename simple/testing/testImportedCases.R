## TRY TO SEE HOW THE METHOD BEHAVES ##

############################################
## DATA SIMULATION ##

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

## w <- c(0,.1,.2,.5,2,.5,.2,.1)
## full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4, mu.transv=0.2e-4)
## dat <- full[1:20]
## collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
##save(w, full, dat, collecDates, file="Robjects/data4.RData")

load("Robjects/data5.RData")
plot(dat, main="Data")
############################################

BURNIN <- 2e4

w <- c(0,.1,.2,.5,2,.5,.2,.1)
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4, mu.transv=0.2e-4)
dat <- full[1:30]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")



## TRY SPECIFYING THE RIGHT OUTLIERS, SEE IF IT WORKS ##
initree <- rep(1, dat$n)
initree[is.na(dat$ances)] <- 0
moveances <- !is.na(dat$ances)

system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree=initree, n.iter=2e5, move.ances=moveances))




############################################
## ESTIMATE EVERYTHING - PARALLEL VERSION ##
## run outbreaker
system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=2e5))

system.time(res <- outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=2e5))

## check results ##
plot.chains(res)

## check ancestries
x <- get.TTree.simple(res, burn=5e4)
temp <- x
temp$ances[is.na(temp$ances)] <- 0
temp2 <- dat$ances
temp2[is.na(temp2)] <- 0
temp2[1] <- NA
mean(temp$ances==temp2,na.rm=TRUE)

v.col <- rep("lightblue",length(x$ances))
notOk <- which(temp$ances!=temp2)
if(length(notOk)>0) v.col[notOk] <- "red"
par(mfrow=c(1,1))
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col)


## check frequency of external infections
alpha <- res$chains[,grep("alpha", names(res$chains))]
barplot(apply(alpha,2, function(e) mean(e==0)), main="Proportion of inferred external case")


## check mutation rates
par(mfrow=c(2,1))
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=1e-4, col="blue")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=1e-4, col="blue")



## check kappa
v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$n.gen!=1)
if(length(notOk)>0) v.col[notOk] <- "red"
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction (red=wrong kappa)", vertex.color=v.col, annot="n.gen")



## check Tinf
toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
boxplot(Tinf, col="grey")
points(dat$dates, col="red", pch="x",cex=2)



## check Pi
plot.chains(res, "pi", omit=BURNIN, type="de")
abline(v=1)



## check Phi
plot.chains(res, "phi", omit=BURNIN, type="de")
abline(v=1/20)

############################################
