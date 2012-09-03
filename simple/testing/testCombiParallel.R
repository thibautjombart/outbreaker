## TRY TO SEE HOW THE METHOD BEHAVES ##

############################################
## DATA SIMULATION ##

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

## w <- c(0,.1,.2,.5,2,.5,.2,.1)
## ####w <- c(0,1,1,1,.5,.2,.1)
## full <- simOutbreak(R0=2, infec.curve=w, mu.transi=2e-4, mu.transv=1e-4)
## dat <- full[1:20]
## collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
## save(w, full, dat, collecDates, file="Robjects/data4.RData")

load("Robjects/data4.RData")
plot(dat, main="Data")
############################################







############################################
## ESTIMATE EVERYTHING - PARALLEL VERSION ##
## run outbreaker
system.time(res <- outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="seqTrack", n.iter=5e5))

## check results ##
plot.chains(res)

par(mfrow=c(2,1))
plot.chains(res, "mu1",type="dens", omit=1e5)
abline(v=2e-4, col="blue")
plot.chains(res, "mu2",type="dens", omit=1e5)
abline(v=1e-4, col="blue")

## check ancestries
x <- get.TTree.simple(res)
mean(x$ances==dat$ances,na.rm=TRUE)

v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$ances!=dat$ances)
if(length(notOk)>0) v.col[notOk] <- "red"
par(mfrow=c(1,1))
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col)

## check frequency of external infections
alpha <- res$chains[,grep("alpha", names(res$chains))]
barplot(apply(alpha,2, function(e) mean(e==0)), main="Proportion of inferred external case")




## check kappa
v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$n.gen!=1)
if(length(notOk)>0) v.col[notOk] <- "red"
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction (red=wrong kappa)", vertex.color=v.col, annot="n.gen")



## check Tinf
toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>1e5, toKeep]
boxplot(Tinf, col="grey")
points(dat$dates, col="red", pch="x",cex=2)



## check Pi
plot.chains(res, "pi", omit=1e5, type="de")
abline(v=1)



## check Phi
plot.chains(res, "phi", omit=1e5, type="de")
abline(v=1/20)

############################################
