## SOURCE LIBRARIES ##
library(ggplot2)


## GET POOLED RESULTS ##
source("/home/thibaut/dev/outbreaker/outbreaker-code/simple/simulations/poolResults.R")
x <- poolResults()


## analysis of high mu results ##
x.highmu <- x[x$type=="mu1-0.0002",]

## link outbreak size / tree quality
cor.test(x.highmu$prop.ances.ok, x.highmu$n)

## link mutation rate / tree quality
muOK <- x.highmu$mu1.ok & x.highmu$mu2.ok
summary(lm(x.highmu$prop.ances.ok~muOK))

## bimodal likelihoods
x.highmu <- x.highmu[order(x.highmu$prop.ances.ok),]

## function to get the likelihood distribution of a simulation
keys <- sub(".$","",rownames(x.highmu))

f1 <- function(key,type="series",what="post"){
    load(paste("mu1-0.0002/",key,"/",key,".RData",sep=""))
    plot.chains(res, what=what, main="likelihood", type=type)
    mtext(side=3, text=paste("simulation:",key,sep=""))
}


## SOME FIGURES ##
## nb of each simulation
table(x$type)

## general results - consensus ancestry OK
qplot(prop.ances.ok, data=x, color=type, geom="density")

## general results - actual ancestor support
qplot(msup.ances, data=x, color=type, geom="density")

## support for actual kappa
qplot(msup.kappa, data=x, color=type, geom="density")

## support for mu1
qplot(mu1.ok, data=x, fill=type, geom="bar")

## support for mu2
qplot(mu2.ok, data=x, fill=type, geom="bar")

## error about the date
qplot(merr.date, data=x, color=type, geom="density")

## support for mu2
qplot(prop.imp.ok, data=x, fill=type, geom="bar")

## imported cases
qplot(prop.imp.ok, data=x, fill=type, geom="bar")

## time vs n
qplot(n, time/60, data=x, color=type, geom=c("point"), ylab="time (minutes)")
