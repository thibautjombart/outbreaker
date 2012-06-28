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
wdens <- c(0.5,1,2,1,0)

## shit density function
##wdens <- c(1,1,1,1,1)

## run outbreaker
res <- outbreaker(as.DNAbin(gen), dates, wdens, n.iter=1e6, sample.every=1000, init.tree="star",quiet=TRUE)


## read outputs
chains <- res$chains
param <- read.table("mcmcOutput.txt",header=TRUE)

## check mu1/mu2
matplot(chains$step, chains[,2:3], type="l", col=c("blue","red"), lty=1, ylim=c(0,0.003))
abline(h=c(0.001,0.002), col=c("blue","red"), lty=2)
legend("topright",fill=c("blue","red"),leg=c("transition","transversion"), title="mutation rates")

plot(density(chains$mu1[-(1:100)]),main="mu1 posterior")
plot(density(chains$mu2[-(1:100)]),main="mu2 posterior")

plot(param$step, param$p_accept_mu1, type="l", main="p accept mu",ylim=c(0,1))


## check gamma
plot(chains$step,chains$gamma,type="l", main="gamma")
plot(density(tail(chains$gamma,round(nrow(chains)/2))),main="gamma")

plot(param$step, param$p_accept_gamma, type="l", main="p accept gamma",ylim=c(0,1))


## check kappa
library(UsingR)
kappa <- tail(chains[,18:23],round(nrow(chains)/2))
apply(kappa,2,mean)
violinplot(kappa)
title("kappa posterior")


## check alpha
f1 <- function(x) {table(x)/length(x)}
alpha <- tail(chains[,12:17],round(nrow(chains)/2))
apply(alpha,2,f1)

assign <- sapply(apply(alpha,2,f1),function(e) names(e)[which.max(e)])
assign # this should be 3-1-4-0-3-5

## check infection dates
infdates <- tail(chains[,6:11],round(nrow(chains)/2))
violinplot(infdates,ylim=c(-3,5),ylab="date")
title("Infection dates posterior")
points(1:6, c(2,4,0,-2,2,4), col="red",lwd=3,cex=1)
mtext(3, text="red = actual date")

plot(param$step, param$p_accept_Tinf, type="l",main="acceptance Tinf",ylim=c(0,1))


## check mcmc param
head(param,50)
tail(param,10)


## plot results ##
x <- get.TTree.simple(res)
x
plot(x, main="Reconstructed outbreak")

}