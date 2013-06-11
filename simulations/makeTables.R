## SOURCE LIBRARIES ##
library(ggplot2)


## GET POOLED RESULTS ##
source("/home/thibaut/dev/outbreaker/outbreaker-code/simple/simulations/poolResults.R")
x <- poolResults()


## SPLIT RESULTS BETWEEN BASIC AND FIXED-MU RESULTS ##
x.basic <- x[-grep("fixedMu", x$type),]
x.fixedMu <- x[grep("fixedMu", x$type),]


## get ordered levels for types of simulation
lev.basic <- sort(unique(gsub("-fixedMu","",x.basic$type)))
type.basic <- factor(x.basic$type, levels=lev.basic) # new factor
all(as.character(type.basic)==x.basic$type) # sanity check - must be TRUE

x.fixedMu$type <- gsub("-fixedMu","",x.fixedMu$type)
lev.fixedMu <- sort(unique(x.fixedMu$type))
type.fixedMu <- factor(x.fixedMu$type, levels=lev.fixedMu) # new factor
all(as.character(type.fixedMu)==x.fixedMu$type) # sanity check - must be TRUE



## MAKE TABLE OF RESULTS ##
f1 <- function(vec,digit=2){
    stat <- c(mean(vec,na.rm=TRUE),quantile(vec, c(.05,.95),na.rm=TRUE))
    stat <- round(stat,digit)
    out <- paste(stat[1], "\n[",stat[2],";",stat[3],"]",sep="")
    return(out)
}


## GET ONLY RELEVANT (REPORTED) RESULTS ##
## BASIC RESULTS ##
res.basic <- x.basic[,c("prop.ances.ok", "msup.ances", "msup.kappa", "merr.date", "prop.imp.ok", "mu1.ok", "mu2.ok", "merr.pi" )]

## get table
tabres.basic <- apply(res.basic, 2, tapply, type.basic, f1)

## rename columns
colnames(tabres.basic) <- c("proportion of correct ancestry","mean support for true ancestor","mean support for true kappa", "mean infection date error", "proportion of imported case detected", "proportion of correct mu1", "proportion of correct mu2", "mean error in pi")

## save table
write.csv(tabres.basic, file="tables/tabRes.basic.csv")




## FIXEDMU RESULTS ##
res.fixedMu <- x.fixedMu[,c("prop.ances.ok", "msup.ances", "msup.kappa", "merr.date", "prop.imp.ok", "merr.pi" )]

## get table
tabres.fixedMu <- apply(res.fixedMu, 2, tapply, type.fixedMu, f1)

## rename columns
colnames(tabres.fixedMu) <- c("proportion of correct ancestry","mean support for true ancestor","mean support for true kappa", "mean infection date error", "proportion of imported case detected", "mean error in pi")

## save table
write.csv(tabres.fixedMu, file="tables/tabRes.fixedMu.csv")

