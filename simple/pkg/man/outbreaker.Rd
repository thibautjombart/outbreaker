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

## dates (collection; infection = 2 days before collection)
dates <- c(4,6,2,0,4,6)

## density function
wdens <- c(0.5,0.5,2,0.2,0)

## shit density function
##wdens <- c(1,1,1,1,1)

## run outbreaker
outbreaker(as.DNAbin(gen), dates, wdens, n.iter=5e5, sample.every=1000, init.tree="star",quiet=TRUE)


## read chains outputs
res <- read.table("output.txt",header=TRUE)
param <- read.table("mcmcOutput.txt",header=TRUE)

## check mu1/mu2
matplot(res$step, res[,2:3], type="l", col=c("blue","red"), lty=1, ylim=c(0,0.003))
abline(h=c(0.001,0.002), col=c("blue","red"), lty=2)
legend("topright",fill=c("blue","red"),leg=c("transition","transversion"), title="mutation rates")

plot(density(res$mu1[-(1:100)]),main="mu1 posterior")
plot(density(res$mu2[-(1:100)]),main="mu2 posterior")

matplot(param$step, param[,3], type="l", main="p accept mu")


## check gamma
plot(res$step,res$gamma,type="l", main="gamma")
plot(density(tail(res$gamma,round(nrow(res)/2))),main="gamma")

matplot(param$step, param[,4], type="l", main="p accept gamma")


## check kappa
library(UsingR)
kappa <- tail(res[,17:22],round(nrow(res)/2))
apply(kappa,2,mean)
violinplot(kappa)
title("kappa")


## check alpha
f1 <- function(x) {table(x)/length(x)}
alpha <- tail(res[,11:16],round(nrow(res)/2))
apply(alpha,2,f1)

assign <- sapply(apply(alpha,2,f1),function(e) names(e)[which.max(e)])
assign # this should be 3-1-4-0-3-5

## check infection dates
infdates <- tail(res[,5:10],round(nrow(res)/2))
violinplot(infdates)
title("Infection dates")
points(1:6, c(2,4,0,-2,2,4), col="red",lwd=3,cex=1)


## check mcmc param
head(param,50)
tail(param,10)


}