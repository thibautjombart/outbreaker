\encoding{UTF-8}
\name{outbreaker}
\alias{outbreaker}
\title{Outbreaker}
\description{
  These functions are under development. Please contact the author if
  you would like to use them.

  \code{outbreaker}: basic implementation of the model.
  \code{outbreaker.parallel}: parallelized version, requiring the
  package \code{parallel}.
}
\usage{
outbreaker(dna, dates, w.dens, w.trunc=length(w.dens),
           init.tree=c("seqTrack","random","star","none"),
           n.iter=1e5, sample.every=500, tune.every=500,
           burnin=2e4, find.import=TRUE, find.import.n=50,
           pi.param1=10, pi.param2=1,
           init.mu1=1e-5, init.gamma=1,
           move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
           move.Tinf=TRUE, move.pi=TRUE,
           quiet=TRUE, res.file.name="output.txt",
           tune.file.name="mcmcOutput.txt", seed=NULL)


outbreaker.parallel(n.runs, parallel=require("parallel"), n.cores=NULL,
                   dna, dates, w.dens, w.trunc=length(w.dens),
                   init.tree=c("seqTrack","random","star","none"),
                   n.iter=1e5, sample.every=500, tune.every=500,
                   burnin=2e4, find.import=TRUE, find.import.n=50,
                   pi.param1=10, pi.param2=1,
                   init.mu1=1e-5, init.gamma=1,
                   move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                   move.Tinf=TRUE, move.pi=TRUE,
                   quiet=TRUE, res.file.name="output.txt",
                   tune.file.name="mcmcOutput.txt", seed=NULL)
}
\arguments{
  \item{dna}{the DNA sequences in \code{DNAbin} format (see
    \code{read.dna} in the ape package}
  \item{dates}{a vector indicating the collection dates, provided either
    as integer numbers or in POSIXct format. Higher numbers should
    correspond to later dates.}
  \item{w.dens}{a vector of numeric values indicating the generation
    time distribution, which indicate how infectious a
    patient is at t=0, 1, 2, ... time steps after infection. By
    default, we impose w.dens[1]=0, meaning that an infected patient
    cannot be infectious at the same time step. If not standardized,
    this distribution is rescaled to sum to 1.}
  \item{w.trunc}{an integer indicating where the distribution
    \code{w.dens} should be truncated, so that further
    values are zero. Length of \code{w.dens} by default.}
  \item{init.tree}{the tree used to initialize the MCMC. Can be either a
    character string indicating how this tree should be computed, or a
    vector if integers corresponding to the tree itself, where the i-th value
    correspond to the index of the ancestor of 'i' (i.e., \code{init.tree[i]}
    is the ancestor of \code{i}. See details.}
  \item{n.iter}{an integer indicating the number of chains in the MCMC,
    including the burnin period; defaults to \code{100,000}. This number
    should be subsequently higher than the burnin.}
  \item{sample.every}{an integer indicating the frequency at which to
    output MCMC results, defaulting to 500 (i.e., output to file every
    500 iterations).}
  \item{tune.every}{an integer indicating the frequency at which
    proposal distributions are tuned, defaulting to 500 (i.e., tune
    proposal distribution every 500 iterations).}
  \item{burnin}{an integer indicating the number of iterations of the
    burnin, after which the chains are supposed to have mixed; parameter
    values are only relevant after the burnin period. Used only when
    imported cases are automatically detected.}
  \item{find.import}{a logical indicating whether imported cases should
    be automatically detected; default to \code{TRUE}. As the model
    assumes that all imported cases are known, this should be left to
    \code{TRUE} unless one is sure that only the first case was imported.}
  \item{find.import.n}{an integer indicating whether how many chains
    should be used to determine imported cases; note that this
    corresponds to chains that are output after the burnin, so that a
    total of (burnin + output.every*find.import.n) chains will be used
    in the prior step to determine imported cases. Defaults to \code{50}.}
  \item{pi.param1, pi.param1}{two numeric values being the parameters of
    the Beta distribution used as a prior for \eqn{\kappa_i}. This prior
    is Beta(10,1) by default, indicating that a majority of cases are
    likely to have been observed.}
  \item{init.mu1,init.gamma}{initial values for the mutation rates (mu1:
    transitions; gamma x mu1: transversions.}
  \item{move.mut,move.ances,move.kappa,move.Tinf,move.pi}{logicals
    indicating whether the named items should be
    estimated ('moved' in the MCMC), or not, all
    defaulting to TRUE. \code{move.mut} handles both mutation rates,
    while \code{move.ances} specifies whether the tree should be estimated.}
  \item{quiet}{a logical indicating whether messages should be displayed
    to the screen.}
  \item{res.file.name}{a character string indicating the name of the
    file used to store MCMC outputs.}
  \item{tune.file.name}{a character string indicating the name of the
    file used to store MCMC tuning parameters.}
  \item{seed}{an integer used to set the random seed of the C procedures.}
  \item{n.runs}{an integer indicating the number of independent MCMC to
    run, either in parallel (if \code{parallel} is used), or serially
    (otherwise).}
  \item{parallel}{a logical indicating whether the package
    \code{parallel} should be used to run parallelized computations; by
    default, this is used if available. This package is currently
    unavailable on Windows systems.}
  \item{n.cores}{an integer indicating the number of cores to be used
    for parallelized computations; if NULL (default value), then all the
    available cores are used.}
}
\value{
  Both procedures return a list with the following components:
  \itemize{
    \item chains: a data.frame containing MCMC outputs (which are also
    stored in the file indicated in \code{res.file.name}).
    \item collec.dates: (data) the collection dates.
    \item w: (data) the infectivity curve (argument \code{w.dens})
    \item D: a matrix of genetic distances (in number of mutations)
    between all pairs of sequences.
    \item tune.end: an integer indicating at which iteration the
    proposal auto-tuning procedures all stopped. All chains prior to
    this are to be considered as burnin.
    \item find.import: a logical indicating if imported cases were to be
    automatically detected.
    \item burnin: an integer indicating the pre-defined burnin, used
    when detecting imported cases.
    \item find.import.at: an integer indicating at which iteration of
    the preliminary MCMC imported cases were detected.
    \item n.runs: the number of independent runs used.
    \item call: the matched call.
  }
}
\references{
 Jombart, Cori, Didelot, Cauchemez, Ferguson and Fraser (in prep). A
 Bayesian approach for reconstructing disease outbreaks based on genomic
 data.
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
library(outbreaker)
library(adegenet)
library(ape)

BURNIN <- 1e4

w <- c(0,.1,.2,.5,2,.5,.2,.1)
mu1 <- 1e-4
mu2 <- 0.2e-4
barplot(w/sum(w), names=0:7,  main="Generation time distribution")
full <- list(n=0)
while(full$n<30){
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=mu1,mu.transv=mu2)
}


####################################
## EXAMPLE WITH COMPLETE SAMPLING ##
####################################

dat <- full[1:30]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")

res <- outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=5e4)


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


## check mutation rates
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=mu1, col="blue")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=mu2, col="blue")


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








######################################
## EXAMPLE WITH INCOMPLETE SAMPLING ##
######################################

dat <- dat.old <- full[sort(sample(1:60, 30))] # 30 random cases from the first 100
dat$id <- 1:length(dat$id)
dat$ances <- match(dat.old$ances, dat.old$id)

collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")

res <- outbreaker.parallel(n.runs=4, dna=dat$dna,dates=collecDates,w.dens=w, init.tree="star", n.iter=10e4, pi.param1=3,pi.param2=3)

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
abline(v=mean(dat$ngen>1))


## check mutation rates
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=mu1, col="blue")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=mu2, col="blue")


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

}