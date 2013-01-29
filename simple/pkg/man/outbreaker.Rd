\encoding{UTF-8}
\name{outbreaker}
\alias{outbreaker}
\alias{outbreaker.parallel}
\title{Outbreaker}
\description{
  These functions are under development. Please contact the author if
  you would like to use them.

  \code{outbreaker}: basic implementation of the model.
  \code{outbreaker.parallel}: parallelized version, requiring the
  package \code{parallel}.
}
\usage{
outbreaker(dna=NULL, dates, idx.dna=NULL,
           w.dens, w.trunc=length(w.dens),
           f.dens=w.dens, f.trunc=length(f.dens),
           init.tree=c("seqTrack","random","star","none"),
           init.kappa=NULL,
           n.iter=1e5, sample.every=500, tune.every=500,
           burnin=2e4, find.import=TRUE, find.import.n=50,
           pi.param1=10, pi.param2=1,
           init.mu1=NULL, init.gamma=1,
           move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
           move.Tinf=TRUE, move.pi=TRUE,
           outlier.threshold = 1000,
           quiet=TRUE, res.file.name="chains.txt", tune.file.name="tuning.txt", seed=NULL)


outbreaker.parallel(n.runs, parallel=require("parallel"), n.cores=NULL,
                    dna=NULL, dates, idx.dna=NULL, w.dens, w.trunc=length(w.dens),
                    f.dens=w.dens, f.trunc=length(f.dens),
                    init.tree=c("seqTrack","random","star","none"),
                    init.kappa=NULL,
                    n.iter=1e5, sample.every=500, tune.every=500,
                    burnin=2e4, find.import=TRUE, find.import.n=50,
                    pi.param1=10, pi.param2=1,
                    init.mu1=NULL, init.gamma=1,
                    move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                    move.Tinf=TRUE, move.pi=TRUE, outlier.threshold = 1000,
                    quiet=TRUE, res.file.name="chains.txt", tune.file.name="tuning.txt", seed=NULL)
}
\arguments{
  \item{dna}{the DNA sequences in \code{DNAbin} format (see
    \code{read.dna} in the ape package}
  \item{dates}{a vector indicating the collection dates, provided either
    as integer numbers or in POSIXct format. By convention, zero indicates
    the oldest date.}
  \item{idx.dna}{an integer vector indicating to which case each dna
    sequence in \code{dna} corresponds. Not required if each case has a
    sequence, and the order of the sequences matches that of the cases.}
  \item{w.dens}{a vector of numeric values indicating the generation
    time distribution, reflecting how infectious a patient is at t=0, 1,
    2, ... time steps after infection. By convention, w.dens[1]=0,
    meaning that an newly infected patient cannot be instantaneously
    infectious. If not standardized, this distribution is rescaled to
    sum to 1.}
  \item{w.trunc}{an integer indicating after which time step the distribution
    \code{w.dens} is zero everywhere (i.e., truncated to zero).}
  \item{f.dens}{similar to \code{w.dens}, except that this is the
    distribution of time interval between infection and sampling.}
  \item{f.trunc}{an integer indicating after which time step the distribution
    \code{f.dens} is zero everywhere (i.e., truncated to zero).}
  \item{init.tree}{the tree used to initialize the MCMC. Can be either a
    character string indicating how this tree should be computed, or a
    vector if integers corresponding to the tree itself, where the i-th value
    correspond to the index of the ancestor of 'i' (i.e., \code{init.tree[i]}
    is the ancestor of \code{i}). Accepted character strings are
    "seqTrack" (uses seqTrack output as initialize tree), "random"
    (ancestor randomly selected from preceding cases), "star" (all cases
    coalesce to the first case), and "none" (no ancestry to start with).}
  \item{init.kappa}{as \code{init.tree}, but values indicate the number
    of generations between each case and its most recent
    sampled ancestor.}
  \item{n.iter}{an integer indicating the number of iterations in the MCMC,
    including the burnin period; defaults to \code{100,000}.}
  \item{sample.every}{an integer indicating the frequency at which to
    sample from the MCMC, defaulting to 500 (i.e., output to file every
    500 iterations).}
  \item{tune.every}{an integer indicating the frequency at which
    proposal distributions are tuned, defaulting to 500 (i.e., tune
    proposal distribution every 500 iterations).}
  \item{burnin}{an integer indicating the number of iterations of the
    burnin, after which the chains are supposed to have mixed; parameter
    values are only relevant after the burnin period. Used only when
    imported cases are automatically detected.}
  \item{find.import}{a logical indicating whether imported cases should
    be automatically detected; default to \code{TRUE}.}
  \item{find.import.n}{an integer indicating how many chains
    should be used to determine imported cases; note that this
    corresponds to chains that are output after the burnin, so that a
    total of (burnin + output.every*find.import.n) chains will be used
    in the prior step to determine imported cases. Defaults to \code{50}.}
  \item{pi.param1, pi.param2}{two numeric values being the parameters of
    the Beta distribution used as a prior for \eqn{\kappa_i}. This prior
    is Beta(10,1) by default, indicating that a majority of cases are
    likely to have been observed. Use Beta(1,1) for a flat prior.}
  \item{init.mu1,init.gamma}{initial values for the mutation rates (mu1:
    transitions; gamma x mu1: transversions.}
  \item{move.mut,move.ances,move.kappa,move.Tinf,move.pi}{logicals
    indicating whether the named items should be estimated ('moved' in
    the MCMC), or not, all defaulting to TRUE. \code{move.mut} handles
    both mutation rates, while \code{move.ances} is a vector of length
    'n' indicating which ancestries should be estimated.}
  \item{outlier.threshold}{a numeric value indicating the threshold for
    detecting low likelihood values corresponding to imported
    cases. Outliers have a likelihood \code{outlier.threshold} smaller
    than the average.}
  \item{quiet}{a logical indicating whether messages should be displayed
    to the screen.}
  \item{res.file.name}{a character string indicating the name of the
    file used to store MCMC outputs.}
  \item{tune.file.name}{a character string indicating the name of the
    file used to store MCMC tuning parameters.}
  \item{seed}{an integer used to set the random seed of the C procedures.}
  \item{n.runs}{an integer indicating the number of independent chains to
    run, either in parallel (if \code{parallel} is used), or serially
    (otherwise).}
  \item{parallel}{a logical indicating whether the package
    \code{parallel} should be used to run parallelized computations; by
    default, it is used if available.}
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
    \item w: (data) the generation time distribution (argument \code{w.dens})
    \item f: (data) the distribution of the time to collection (argument \code{f.dens})
    \item D: a matrix of genetic distances (in number of mutations)
    between all pairs of sequences.
    \item idx.dna: (data) the index of the case each dna sequence corresponds to
    \item tune.end: an integer indicating at which iteration the
    proposal auto-tuning procedures all stopped.
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
 Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson N (in prep).
Bayesian reconstruction of disease outbreaks by combining epidemiologic and genomic data.
}
\author{
  See 'references' for the model.
  Implementation by Thibaut Jombart (\email{t.jombart@imperial.ac.uk})
}
\examples{
\dontrun{
library(outbreaker)

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

dat <- full[1:10]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")


#### TEST DIFFERENT VERSIONS OF THE METHOD ####
## test basic ##
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w,init.tree="star", n.iter=3e4, find.import=FALSE)

## test without genetic info ##
res <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, init.tree="star", n.iter=3e4, find.import=FALSE)

## test without import detection ##
res <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, init.tree="star", n.iter=3e4)

## test without f different from w ##
res <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, f.dens=1, init.tree="star", n.iter=3e4)

## test parallel ##
res <- outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=3e4, find.import=FALSE)


#### TEST POST-PROCESSING OF THE RESULTS ####
## test plots ##
plotChains(res)
plot(get.tTree(res))
transGraph(res, thres=0.1)
plotOutbreak(res)

## test incidence curves ##
get.incid(res, main="Incidence curves")

## test effective reproduction numbers ##
barplot(table(get.R(res)), main="Effective reproduction number - distribution")

get.Rt(res, main="Effective reproduction number over time")

}
}