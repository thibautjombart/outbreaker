\encoding{UTF-8}
\name{simulated outbreak dataset}
\alias{toyOutbreak}
\docType{data}
\title{Toy outbreak dataset used to illustrate outbreaker}
\description{
  This toy outbreak dataset was simulated using
  \code{\link{simOutbreak}}. This dataset is a list containing the
  following components:
  \itemize{
    \item dat: the data, output of simOutbreak; see \code{dat$call} for
    the actual command line that was used.
    \item w: the generation time distribution, for t=0 to t=3.
    \item collecDates: simulated collection dates dates.
    \item res: the results of \code{outbreaker.parallel}; see
    \code{res$call} for the actual command line that was used.
  }
}
\usage{data(toyOutbreak)}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\examples{
\dontrun{
## command line to generate similar data
w <- c(0, 0.5, 1, 0.75)
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100)[1:30] #works only if outbreak has at least 30 cases
collecDates <- dat$dates + sample(0:3, size=30, replace=TRUE, prob=w)

## command line to analyse the data (4 runs in parallel)
res <-  outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4)
}

## load data
data(toyOutbreak)
attach(toyOutbreak)

## visualize dynamics
matplot(dat$dynam, type="o", pch=20, lty=1, xlim=c(0,32), main="Outbreak dynamics")
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)

## visualize transmission tree
plot(dat, annot="dist", main="Data - transmission tree")

## visualize NJ tree
if(require(ape)){
  plot(nj(dist.dna(dat$dna)))
  title("NJ tree")
}

## assess convergence of chains
plotChains(res)

## represent posterior ancestries
transGraph(res, annot="", main="Posterior ancestries")

## get consensus ancestries
tre <- get.tTree(res)
plot(tre, annot="", main="Consensus ancestries")

## show discrepancies
col <- rep("lightgrey", 30)
col[which(dat$ances != tre$ances)] <- "pink"
plot(tre, annot="", vertex.color=col, main="Consensus ancestries")
mtext(side=3, text="cases with erroneous ancestries in pink")

detach(toyOutbreak)

}
\keyword{datasets}
