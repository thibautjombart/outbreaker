\encoding{UTF-8}
\name{outbreaker}
\alias{outbreaker}
\alias{outbreaker.parallel}
\title{Outbreaker: disease outbreak reconstruction using genetic data}
\description{
  \code{outbreaker} is a tool for the reconstruction of disease
  outbreaks using pathogens genome sequences. It relies on a
  probabilistic model of disease transmission which takes the genetic
  diversity, collection dates, duration of pathogen colonization and
  time interval between cases into account. It is embedded in a Bayesian
  framework which allows to estimate the distributions of parameters of
  interest. It currently allows to estimate:
  \itemize{
    \item transmission trees
    \item dates of infection
    \item missing cases in a chain of transmission
    \item imported cases
    \item mutation rates
  }

  The function \code{outbreaker} is the basic implementation of the
  model. \code{outbreaker.parallel} allows to run several independent
  MCMC in parallel across different cores / processors of the same
  computer. This requires the base package \code{parallel}.

  The spatial module implemented in outbreaker is currently under
  development.  Please contact the author before using it.

}
\usage{
outbreaker(dna=NULL, dates, idx.dna=NULL, mut.model=1, w.dens,
           w.trunc=length(w.dens), f.dens=w.dens,
           f.trunc=length(f.dens), dist.mat=NULL, spa.model=1,
           init.tree=c("seqTrack","random","star"), init.kappa=NULL,
           init.mu1=NULL, init.mu2=init.mu1, init.spa1=NULL,
           init.spa2=NULL, n.iter=1e5, sample.every=500, tune.every=500,
           burnin=2e4, find.import=TRUE, find.import.n=50, pi.param1=10,
           pi.param2=1, spa1.prior=1, spa2.prior=1, move.mut=TRUE,
           move.ances=TRUE, move.kappa=TRUE, move.Tinf=TRUE,
           move.pi=TRUE, move.spa=TRUE, outlier.threshold = 5,
           quiet=TRUE, res.file.name="chains.txt",
           tune.file.name="tuning.txt", seed=NULL )


outbreaker.parallel(n.runs, parallel=require("parallel"), n.cores=NULL,
                    dna=NULL, dates, idx.dna=NULL, mut.model=1, w.dens,
                    w.trunc=length(w.dens), f.dens=w.dens,
                    f.trunc=length(f.dens), dist.mat=NULL, spa.model=1,
                    init.tree=c("seqTrack","random","star"),
                    init.kappa=NULL, init.mu1=NULL, init.mu2=init.mu1,
                    init.spa1=NULL, init.spa2=NULL, n.iter=1e5,
                    sample.every=500, tune.every=500, burnin=2e4,
                    find.import=TRUE, find.import.n=50, pi.param1=10,
                    pi.param2=1, spa1.prior=1, spa2.prior=1,
                    move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                    move.Tinf=TRUE, move.pi=TRUE, move.spa=TRUE,
                    outlier.threshold = 5, quiet=TRUE,
                    res.file.name="chains.txt",
                    tune.file.name="tuning.txt", seed=NULL)
}
\arguments{
  \item{dna}{the DNA sequences in \code{DNAbin} format (see
    \code{\link[ape]{read.dna}} in the ape package); this can be
    imported from a fasta file (extension .fa, .fas, or .fasta) using
    \code{adegenet}'s function \code{\link[adegenet]{fasta2DNAbin}}.}

  \item{dates}{a vector indicating the collection dates, provided either
    as integer numbers or in a usual date format such as \code{Date} or
    \code{POSIXct} format. By convention, zero will indicate the oldest
    date.}
  
  \item{idx.dna}{an optional integer vector indicating to which case
    each dna sequence in \code{dna} corresponds. Not required if each
    case has a sequence, and the order of the sequences matches that of
    the cases.}

  \item{mut.model}{an integer indicating the mutational model to be
    used; 1: one single mutation rate; 2: two rates,
    transitions/transversions.}
  
  \item{spa.model}{an integer indicating the spatial model to be
    used. 0: no spatial model. 1: exponential kernel.}
  
  \item{w.dens}{a vector of numeric values indicating the generation
    time distribution, reflecting the infectious potential of a case
    t=0, 1, 2, ... time steps after infection. By convention,
    w.dens[1]=0, meaning that an newly infected patient cannot be
    instantaneously infectious. If not standardized, this distribution
    is rescaled to sum to 1.}
  
  \item{w.trunc}{an integer indicating after which time step the
    distribution \code{w.dens} should be truncated to zero (effectively,
    the maximum duration of the infectious period).}
  
  \item{f.dens}{similar to \code{w.dens}, except that this is the
    distribution of the colonization time, i.e. time interval during
    which the pathogen can be sampled from the patient.}
  
  \item{f.trunc}{an integer indicating after which time step the
    distribution \code{f.dens} should be truncated to zero (effectively,
    the maximum duration of colonization).}

  \item{dist.mat}{a matrix of pairwise spatial distances between the
    cases.}
  
  \item{init.tree}{the tree used to initialize the MCMC. Can be either a
    character string indicating how this tree should be computed, or a
    vector if integers corresponding to the tree itself, where the i-th
    value corresponds to the index of the ancestor of 'i' (i.e.,
    \code{init.tree[i]} is the ancestor of case \code{i}). Accepted
    character strings are "seqTrack" (uses seqTrack output as initialize
    tree), "random" (ancestor randomly selected from preceding cases),
    and "star" (all cases coalesce to the first case). Note that for
    SeqTrack, all cases should have been sequenced.}
  
  \item{init.kappa}{as \code{init.tree}, but values indicate the number
    of generations between each case and its most recent sampled
    ancestor.}

  \item{n.iter}{an integer indicating the number of iterations in the
    MCMC, including the burnin period; defaults to \code{100,000}.}
  
  \item{sample.every}{an integer indicating the frequency at which to
    sample from the MCMC, defaulting to 500 (i.e., output to file every
    500 iterations).}
  
  \item{tune.every}{an integer indicating the frequency at which
    proposal distributions are tuned, defaulting to 500 (i.e., tune
    proposal distribution every 500 iterations).}
  
  \item{burnin}{an integer indicating the number of iterations for the
    burnin period, after which the chains are supposed to have mixed;
    estimated values of parameter are only relevant after the burnin
    period. Used only when imported cases are automatically detected.}
  
  \item{find.import}{a logical indicating whether imported cases should
    be automatically detected; defaults to \code{TRUE}.}
  
  \item{find.import.n}{an integer indicating how many chains should be
    used to determine imported cases; note that this corresponds to
    chains that are output after the burnin, so that a total of (burnin
    + output.every*find.import.n) chains will be used in the prior run
    to determine imported cases. Defaults to \code{50}.}
  
  \item{pi.param1, pi.param2}{two numeric values being the parameters of
    the Beta distribution used as a prior for \eqn{\kappa_i}. This prior
    is Beta(10,1) by default, indicating that a majority of cases are
    likely to have been observed. Use Beta(1,1) for a flat prior.}
  
  \item{init.mu1,init.mu2}{initial values for the mutation rates (mu1:
    transitions; mu2: transversions).}
  

  \item{init.spa1,init.spa2}{initial values of the spatial parameters.}

  \item{spa1.prior,spa2.prior}{parameters of the prior distribution for
    the spatial parameters. In the spatial model 1, \code{spa1.prior} is
    the mean of an exponential distribution.}

  \item{move.mut,move.ances,move.kappa,move.Tinf,move.pi,move.spa}{logicals
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
  
  \item{seed}{an integer used to set the random seed of the C
    procedures.}  \item{n.runs}{an integer indicating the number of
    independent chains to run, either in parallel (if \code{parallel} is
    used), or serially (otherwise).}
  
  \item{parallel}{a logical indicating whether the package
    \code{parallel} should be used to run parallelized computations; by
    default, it is used if available.}
  
  \item{n.cores}{an integer indicating the number of cores to be used
    for parallelized computations; if NULL (default value), then up to 6
    cores are used, depending on availability.}
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
 Jombart T, Cori A, Didelot X, Cauchemez S, Fraser C and Ferguson N
(accepted).  Bayesian reconstruction of disease outbreaks by combining
epidemiologic and genomic data. PLoS Computational Biology.

}
\author{
  Thibaut Jombart (\email{t.jombart@imperial.ac.uk})
}
\examples{
\dontrun{

#### TEST DIFFERENT VERSIONS OF THE METHOD ####
## test basic ##
res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, n.iter=5e4)
plotChains(res)
par(mfrow=c(1,2))
plot(dat)
title("Simulated data")
plot(get.tTree(res))
title("Reconstructed outbreak")

## test without genetic info ##
res <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, n.iter=5e4, find.import=FALSE)

## test without import detection ##
res <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, init.tree="star", n.iter=3e4)

## test without f different from w ##
res <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, f.dens=1,
init.tree="star", n.iter=3e4)

## test parallel ##
res <- outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates,
w.dens=w, init.tree="star", n.iter=3e4, find.import=FALSE)



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