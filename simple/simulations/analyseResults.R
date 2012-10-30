## SOURCE LIBRARIES ##
library(ggplot2)


## GET POOLED RESULTS ##
source("/home/thibaut/dev/outbreaker/outbreaker-code/simple/simulations/poolResults.R")
x <- poolResults()


## ANALYSE RESULTS ##
## general results - consensus ancestry OK
qplot(prop.ances.ok, data=x, color=type, geom="density")


## general results - actual ancestor support
qplot(msup.ances, data=x, color=type, geom="density")


## support for actual kappa
qplot(msup.kappa, data=x, color=type, geom="density")

## support for mu1
qplot(mu1.ok, data=x, fill=type, geom="bar")

## support for mu2
qplot(mu2.ok, data=x, fill=type, geom="bar")

## error about the date
qplot(merr.date, data=x, color=type, geom="density")

## support for mu2
qplot(prop.imp.ok, data=x, fill=type, geom="bar")

## imported cases
qplot(prop.imp.ok, data=x, fill=type, geom="bar")

## time vs n
qplot(n, time, data=x, color=type, geom=c("point"))
