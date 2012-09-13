library(outbreaker)
library(adegenet)
library(ape)

BURNIN <- 1e4

w <- c(0,.1,.2,.5,2,.5,.2,.1)
barplot(w/sum(w), names=0:7,  main="Generation time distribution")
full <- list(n=0)
while(full$n<30){
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4,mu.transv=0.2e-4)
}
dat <- full[1:30]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")

## RUN OUTBREAKER ##
system.time(res <- outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="seqTrack", n.iter=6e4))


pdf("figures/expl-posterior.pdf")
plot.chains(res, main="Posterior values")
dev.off()


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

pdf("figures/expl-data.pdf")
par(mfrow=c(1,1))
plot(dat,main="data", vertex.color=v.col)
dev.off()

pdf("figures/expl-recons.pdf")
plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col, annot="prob")
mtext(side=3, text=" - figures on arrows indicate posterior support -")
dev.off()


## check mutation rates
par(xpd=FALSE)
pdf("figures/expl-mu1.pdf")
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=1e-4, col="blue")
dev.off()

pdf("figures/expl-mu2.pdf")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=0.2e-4, col="blue")
dev.off()



## check infection dates
toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
colnames(Tinf) <- sub("Tinf_", "case ",colnames(Tinf))

pdf("figures/expl-Tinf.pdf")
boxplot(Tinf, col=funky(20), horizontal=TRUE, las=1, xlab="Time", main="Inference of infection dates")
mtext(side=3, text="(X = actual date)")
points(dat$dates,1:dat$n, col="black", pch="x",cex=2.5)
points(dat$dates,1:dat$n, col=funky(20), pch="x",cex=2)
dev.off()


## get incidence curves
pdf("figures/expl-incidence.pdf")
get.incid(res, main="Incidence curves")
dev.off()


## get effective reproduction numbers
pdf("figures/expl-Rt.pdf")
get.Rt(res, main="Effective reproduction number")
dev.off()
