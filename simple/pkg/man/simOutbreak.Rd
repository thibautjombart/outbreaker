\name{simOutbreak}
\alias{simOutbreak}
\alias{print.simOutbreak}
\alias{[.simOutbreak}
\alias{labels.simOutbreak}
\alias{simOutbreak-class}
\alias{as.igraph.simOutbreak}
\alias{plot.simOutbreak}
\title{Simulation of pathogen genotypes during disease outbreaks}
\description{
  The function \code{simOutbreak} implements simulations of disease
  outbreaks. This functions are currently under development. Please
  email the author before using them.
}
\usage{
simOutbreak(R0, infec.curve, n.hosts = 200, duration = 50, seq.length = 10000, 
    mu.transi = 1e-04, mu.transv = mu.transi/2, rate.import.case = 0.1, 
    diverg.import = 20)
\method{print}{simOutbreak}(x, \dots)
\method{[}{simOutbreak}(x, i, j, drop=FALSE)
\method{labels}{simOutbreak}(object, \dots)
\method{as.igraph}{simOutbreak}(x, \dots)
\method{plot}{simOutbreak}(x, y=NULL, cex=1, col=num2col(x$dates), label=x$id,
        edge.col=num2col(x$nmut[-1], col.pal=seasun), lwd=1, \dots)
}
\arguments{
  \item{R0}{the basic reproduction number}
  \item{infec.curve}{a \code{numeric} vector describing the
individual infectiousness at time t=0, 1, \dots}
  \item{n.hosts}{the number of susceptible hosts at the begining of the
outbreak}
  \item{duration}{the number of time steps for which simulation is run}
  \item{seq.length}{an integer indicating the length of the simulated
    haplotypes, in number of nucleotides.}
  \item{mu.transi}{the rate of transitions, in number of mutation per site and per
    time unit.}
  \item{mu.transv}{the rate of transversions, in number of mutation per site and per
    time unit.}
%   \item{tree}{a \code{logical} indicating whether a (minimum evolution)
% tree of the outbreak should be included in the output.}
  \item{rate.import.case}{the rate at which cases are imported at each
time step.}
  \item{diverg.import}{the number of time steps to the MRCA of all
imported cases.}
  \item{x,object}{\code{simOutbreak} objects.}
  \item{i,j, drop}{\code{i} is a vector used for subsetting the object. For
    instance, \code{i=1:3} will retain only the first three haplotypes of the
    outbreak. \code{j} and \code{drop} are only provided for compatibility,
    but not used.}
  \item{y}{present for compatibility with the generic 'plot'
method. Currently not used.}
  \item{cex}{a size factor for the vertices of the plotted graph.}
  \item{col}{the color of the vertices of the plotted graph.}
  \item{label}{the labels of the vertices of the plotted graph.}
  \item{edge.col}{the color of the edges of the plotted graph.}
  \item{lwd}{the width of the edges of the plotted graph.}
  \item{\dots}{further arguments to be passed to other methods}

}
\author{
Implementation by Thibaut Jombart \email{t.jombart@imperial.ac.uk}
Epidemiological model designed by Anne Cori.}
\value{
  === simOutbreak class ===\cr
  \code{simOutbreak} objects are lists containing the following
slots:\cr
\itemize{
  \item n: the number of cases in the outbreak\cr
  \item dna: DNA sequences in the DNAbin matrix format\cr
  \item dates: infection dates\cr
  \item dynam: a data.frame containing, for each time step (row), the number
of susceptible, infected, or recovered in the population. \cr
  \item id: a vector of integers identifying the cases\cr
  \item ances: a vector of integers identifying infectors
('ancestor')\cr
  \item nmut: the number of mutations corresponding to each ancestry\cr
  \item ngen: the number of generations corresponding to each ancestry\cr
  \item call: the matched call
}
}
\examples{
%%\dontrun{
if(require(ape)&&require(igraph)){
dat <- simOutbreak(R0 = 2, infec.curve = c(0, 1, 1, 1), n.hosts = 100)
dat

## plot first 30 cases
N <- dat$n
plot(dat[1:(min(N,30))], main="First 30 cases")
mtext(side=3, text="nb mutations / nb generations")

## plot a random subset (n=10) of the first cases
x <- dat[sample(1:min(N,30), 10, replace=FALSE)]
plot(x, main="Random sample of 10 of the first 30 cases")
mtext(side=3, text="nb mutations / nb generations")

## plot population dynamics
head(dat$dynam,15)
matplot(dat$dynam[1:max(dat$dates),],xlab="time",ylab="nb of individuals",pch=c("S","I","R"), type="b")
%%}
}
}