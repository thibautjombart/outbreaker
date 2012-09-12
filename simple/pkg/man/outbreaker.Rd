\encoding{UTF-8}
\name{outbreaker}
\alias{outbreaker}
\title{Outbreaker}
\description{
 ...
}
\usage{
}
\arguments{
  \item{x}{...}
}
\value{
  Something cool.
}

\references{
 Us.
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
library(outbreaker)
library(adegenet)
library(ape)

BURNIN <- 1e4

w <- c(0,.1,.2,.5,2,.5,.2,.1)
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4, mu.transv=0.2e-4)
dat <- full[1:20]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")

res <- outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=4e5)


plot.chains(res)

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
x11();
plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col)


## check mutation rates
par(mfrow=c(2,1))
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=1e-4, col="blue")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=1e-4, col="blue")

## check infection dates
toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
colnames(Tinf) <- sub("Tinf_", "case ",colnames(Tinf))
boxplot(Tinf, col=funky(20), horizontal=TRUE, las=1, xlab="Time", main="Inference of infection dates")
mtext(side=3, text="(X = actual date)")
points(dat$dates,1:20, col="black", pch="x",cex=2.5)
points(dat$dates,1:20, col=funky(20), pch="x",cex=2)


## get incidence curves
get.incid(res)



## get effective reproduction numbers
get.Rt(res)

}