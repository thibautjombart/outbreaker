
library(outbreaker)

## GENERATE DATA ##
BURNIN <- 1e4

w <- c(0,.1,.2,.5,2,.5,.2,.1)
barplot(w/sum(w), names=0:7,  main="Generation time distribution")
full <- list(n=0)
while(full$n<30){
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4,mu.transv=0.2e-4)
}


dat <- dat.old <- full[sort(sample(1:60, 30))] # 30 random cases from the first 100
dat$id <- 1:length(dat$id)
dat$ances <- match(dat.old$ances, dat.old$id)



collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")
mtext(side=3, "# mut / # generations")


## RUN OUTBREAKER ##
## fixing kappa to true values where kappa>1
##
kappaToMove <- dat$ngen < 2
kappaInit <- rep(0,dat$n)
kappaInit[!kappaToMove] <- dat$ngen[!kappaToMove]
alphaInit <- c(0, rep(1,dat$n-1))
alphaInit[!kappaToMove] <- dat$ances[!kappaToMove]


res <- outbreaker.parallel(n.runs=4, dna=dat$dna,dates=collecDates,w.dens=w, init.tree=alphaInit, n.iter=6e4,
                           init.kappa=kappaInit, move.ances=kappaToMove, move.kappa=kappaToMove)

plot.chains(res, main="Posterior probabilities")

## check ancestries
x <- get.TTree.simple(res, burn=BURNIN)
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
plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col)


## check nb of generations
temp <- x
temp$n.gen[is.na(temp$ances)] <- NA
temp2 <- dat$ngen
temp2[is.na(temp$ances)] <- NA
temp2[1] <- NA
mean(temp$n.gen==temp2,na.rm=TRUE)

v.col <- rep("lightblue",length(x$ances))
notOk <- which(temp$n.gen!=temp2)
if(length(notOk)>0) v.col[notOk] <- "red"
par(mfrow=c(1,1))
plot(dat,main="data", vertex.color=v.col)
plot(x,main="reconstruction (red=wrong number of generations)", vertex.color=v.col)


## check pi
plot.chains(res, "pi",type="dens", omit=2e4)
abline(v=mean(dat$ngen==1))


## check mutation rates
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=1e-4, col="blue")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=0.2e-4, col="blue")

## check infection dates
toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
colnames(Tinf) <- sub("Tinf_", "case ",colnames(Tinf))
boxplot(Tinf, col=funky(20), horizontal=TRUE, las=1, xlab="Time", main="Inference of infection dates")
mtext(side=3, text="(X = actual date)")
points(dat$dates,1:dat$n, col="black", pch="x",cex=2.5)
points(dat$dates,1:dat$n, col=funky(20), pch="x",cex=2)


## get incidence curves
get.incid(res, main="Incidence curves")


## get effective reproduction numbers
get.Rt(res, main="Effective reproduction number")
