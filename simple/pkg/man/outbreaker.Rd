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
## generate toy dataset ##
##
##        -> 1 -> 2
##       |
##  4 -> 3
##       |
##        -> 5 -> 6
##
## transi: 0.001  transver: 0.002 (per site and generation)
## generation time ~ Poisson(2)

library(ape)

## genomes
gen <- matrix(rep(rep(c('a','g','c','t'), each=250),6),byrow=TRUE,
nrow=6)

gen[3,c(1,751,752)] <- c('g','a','a')

gen[1,] <- gen[3,]
gen[1,c(2,251,252)] <- c('g','c','c')

gen[2,] <- gen[1,]
gen[2,c(3,253,254)] <- c('g','c','c')

gen[5,] <- gen[3,]
gen[5,c(255,501,502)] <- c('a','a','a')

gen[6,] <- gen[5,]
gen[6,c(256,503,504)] <- c('a','g','g')

##plot(nj(dist.dna(as.DNAbin(gen), model="raw")), type="unr")

## dates
dates <- c(4,6,2,0,4,6)


## run outbreaker
outbreaker(as.DNAbin(gen), dates, n.iter=5e5, sample.every=500, sigma.mu1=1e-3)


## read chains outputs
res <- read.table("output.txt",header=TRUE)

matplot(res$step, res[,2:3], type="l", col=c("blue","red"), lty=1, ylim=c(0,0.003))
abline(h=c(0.001,0.002), col=c("blue","red"), lty=2)
legend("topright",fill=c("blue","red"),leg=c("transition","transversion"), title="mutation rates")

plot(density(res$mu1[-(1:100)]),main="mu1 posterior")
plot(density(res$mu2[-(1:100)]),main="mu2 posterior")


## check mcmc parameters
param <- read.table("mcmcOutput.txt",header=TRUE)
matplot(param$step, param[,3], type="l", main="p accept mu")
head(param,50)
tail(param,10)

}