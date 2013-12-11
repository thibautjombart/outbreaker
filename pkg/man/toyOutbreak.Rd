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
## COMMAND LINES TO GENERATE SIMILAR DATA ##
w <- c(0, 0.5, 1, 0.75)
## note: this works only if outbreak has at least 30 case
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100)[1:30]
collecDates <- dat$onset + sample(0:3, size=30, replace=TRUE, prob=w)

## COMMAND LINE TO ANALYSE THE DATA (4 RUNS IN PARALLEL)
res <-  outbreaker.parallel(n.runs=4, dna=dat$dna,
   dates=collecDates,w.dens=w, n.iter=5e4)



## EXAMPLE USING TOYOUTBREAK ##
## LOAD DATA, SET RANDOM SEED
data(toyOutbreak)
attach(toyOutbreak)
set.seed(1)

## VISUALIZE DYNAMICS
matplot(dat$dynam, type="o", pch=20, lty=1,
   main="Outbreak dynamics")
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)

## VISUALIZE TRANSMISSION TREE
plot(dat, annot="dist", main="Data - transmission tree")
mtext(side=3, "arrow annotations are numbers of mutations")


## RUN OUTBREAKER - PARALLEL VERSION
## (takes < 1 min))
res <-  outbreaker.parallel(n.runs=4, dna=dat$dna,
   dates=collecDates,w.dens=w, n.iter=5e4)

## ASSESS CONVERGENCE OF CHAINS
plotChains(res)

## REPRESENT POSTERIOR ANCESTRIES
transGraph(res, annot="", main="Posterior ancestries", thres=.01)

## GET CONSENSUS ANCESTRIES
tre <- get.tTree(res)
plot(tre, annot="", main="Consensus ancestries")

## SHOW DISCREPANCIES
col <- rep("lightgrey", 30)
col[which(dat$ances != tre$ances)] <- "pink"
plot(tre, annot="", vertex.color=col, main="Consensus ancestries")
mtext(side=3, text="cases with erroneous ancestries in pink")

detach(toyOutbreak)
}
}
\keyword{datasets}
