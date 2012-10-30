## source libraries ##
library(ggplot2)


## get pooled results ##
source("/home/thibaut/dev/outbreaker/outbreaker-code/simple/simulations/poolResults.R")
x <- poolResults()


## basic plot ##
qplot(prop.ances.ok, data=x, color=type, geom="density")

