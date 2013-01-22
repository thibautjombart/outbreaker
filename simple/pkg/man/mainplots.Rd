\name{plot.chains}
\alias{plot.chains}
\title{Plot outbreaker's MCMCs}
\description{
  The function \code{plot.chains} is used for plotting MCMCs output by
  \code{outbreaker} and \code{outbreaker.parallel}.
}
\usage{
plot.chains(x, what="post", type=c("series","density"), omit.first=0,
            dens.all=TRUE, col=rainbow(x$n.runs), lty=1, lwd=1, main=what, \dots)
}
\arguments{
  \item{x}{the output of \code{outbreaker} or \code{outbreaker.parallel}.}
  \item{what}{a character chains giving the name of the item to be
    plotted. See \code{names(x$chains)} for possible values. By default,
    log-posterior values are plotted}
  \item{type}{a character indicating if the chains should be plotted as
    time series ("series"), or as density ("density").}
  \item{omit.first}{an integer indicating the number of MCMC steps to
    discard before plotting chains.}
  \item{dens.all}{a logical indicating if, in the case of multiple runs,
    the overall density of the different chains should be plotted in
    addition to individual densities.}
  \item{col}{a vector of colors to be used to plot different chains.}
  \item{lty}{a vector of integers specifying line types for the
    different chains.}
  \item{lwd}{same as \code{lty}, but for line width.}
  \item{main}{the title to be added to the plot.}
  \item{\dots}{further arguments to be passed to the plot function.}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{

}
