pkgname <- "outbreaker"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('outbreaker')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rstat")
### * Rstat

flush(stderr()); flush(stdout())

### Name: reproduction numbers
### Title: Derive reproduction numbers from outbreak's outputs
### Aliases: get.R get.Rt get.incid

### ** Examples

## load data
data(fakeOutbreak)
attach(fakeOutbreak)

## individual R
barplot(table(get.R(res)), main="Individual effective reproduction numbers")

## R(t)
get.Rt(res)

## incidence
get.incid(res)

detach(fakeOutbreak)




cleanEx()
nameEx("fakeOutbreak")
### * fakeOutbreak

flush(stderr()); flush(stdout())

### Name: simulated outbreak dataset
### Title: Toy outbreak dataset used to illustrate outbreaker
### Aliases: fakeOutbreak
### Keywords: datasets

### ** Examples

## Not run: 
##D ## COMMAND LINES TO GENERATE SIMILAR DATA ##
##D w <- c(0, 0.5, 1, 0.75)
##D ## note: this works only if outbreak has at least 30 case
##D dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100)[1:30]
##D collecDates <- dat$onset + sample(0:3, size=30, replace=TRUE, prob=w)
## End(Not run)

## EXAMPLE USING TOYOUTBREAK ##
## LOAD DATA, SET RANDOM SEED
data(fakeOutbreak)
attach(fakeOutbreak)

## VISUALIZE DYNAMICS
matplot(dat$dynam, type="o", pch=20, lty=1,
   main="Outbreak dynamics", xlim=c(0,28))
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)

## VISUALIZE TRANSMISSION TREE
plot(dat, annot="dist", main="Data - transmission tree")
mtext(side=3, "arrow annotations are numbers of mutations")


## Not run: 
##D ## RUN OUTBREAKER - PARALLEL VERSION
##D ## (takes < 1 min))
##D set.seed(1)
##D res <-  outbreaker.parallel(n.runs=4, dna=dat$dna,
##D    dates=collecDates,w.dens=w, n.iter=5e4)
## End(Not run)


## ASSESS CONVERGENCE OF CHAINS
plotChains(res)
plotChains(res, burnin=2e4)

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

## GET EFFECTIVE REPRODUCTION OVER TIME
get.Rt(res)

## GET INDIVIDUAL EFFECTIVE REPRODUCTION
head(get.R(res))
boxplot(get.R(res), col="grey", xlab="Case",
        ylab="Effective reproduction number")

## GET MUTATION RATE PER TIME UNIT
## per genome
head(get.mu(res))

## per nucleotide
mu <- get.mu(res, genome.size=1e4)
head(mu)

summary(mu)
hist(mu, border="lightgrey", col="grey", xlab="Mutation per day and nucleotide",
     main="Posterior distribution of mutation rate")

detach(fakeOutbreak)




cleanEx()
nameEx("getMu")
### * getMu

flush(stderr()); flush(stdout())

### Name: Mutation rate estimation
### Title: Derive mutation rate estimation from outbreak's outputs
### Aliases: get.mu

### ** Examples

## load data
data(fakeOutbreak)
attach(fakeOutbreak)

mu <- get.mu(res, genome.size=ncol(dat$dna))
hist(mu, col="grey",
     main="Inferred distribution of mu",
     xlab="mutations/site/day")
abline(v=1e-4,lty=2, lwd=4, col="royalblue")
mtext(side=3, "Dashed line = actual value")

detach(fakeOutbreak)



cleanEx()
nameEx("mainplots")
### * mainplots

flush(stderr()); flush(stdout())

### Name: outbreaker graphics
### Title: Plot outbreaker's results
### Aliases: plotChains transGraph plotOutbreak

### ** Examples

data(fakeOutbreak)
attach(fakeOutbreak)

## examine MCMC
plotChains(res)
plotChains(res,type="dens")
plotChains(res,type="dens", what="mu1", burnin=2e4)

## represent posterior ancestries
transGraph(res, annot="", main="Posterior ancestries")
transGraph(res, annot="", main="Posterior ancestries - support > 0.5",
   threshold=0.5)
if(require(adegenet)){
transGraph(res, annot="", main="Posterior ancestries - support > 0.01",
   threshold=0.01, col.pal=spectral)
}
## summary plot
plotOutbreak(res,cex.bubble=0.5, thres.hide=0.5,
   main="Outbreak reconstruction")


detach(fakeOutbreak)




cleanEx()
nameEx("outbreaker")
### * outbreaker

flush(stderr()); flush(stdout())

### Name: outbreaker
### Title: Outbreaker: disease outbreak reconstruction using genetic data
### Aliases: outbreaker outbreaker.parallel

### ** Examples


## EXAMPLE USING TOYOUTBREAK ##
## LOAD DATA, SET RANDOM SEED
data(fakeOutbreak)
attach(fakeOutbreak)

## VISUALIZE DYNAMICS
matplot(dat$dynam, type="o", pch=20, lty=1,
   main="Outbreak dynamics", xlim=c(0,28))
legend("topright", legend=c("S","I","R"), lty=1, col=1:3)

## VISUALIZE TRANSMISSION TREE
plot(dat, annot="dist", main="Data - transmission tree")
mtext(side=3, "arrow annotations are numbers of mutations")


## Not run: 
##D ## RUN OUTBREAKER - PARALLEL VERSION
##D ## (takes < 1 min))
##D set.seed(1)
##D res <-  outbreaker.parallel(n.runs=4, dna=dat$dna,
##D    dates=collecDates,w.dens=w, n.iter=5e4)
## End(Not run)


## ASSESS CONVERGENCE OF CHAINS
plotChains(res)
plotChains(res, burnin=2e4)

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

## GET EFFECTIVE REPRODUCTION OVER TIME
get.Rt(res)

## GET INDIVIDUAL EFFECTIVE REPRODUCTION
head(get.R(res))
boxplot(get.R(res), col="grey", xlab="Case",
        ylab="Effective reproduction number")

## GET MUTATION RATE PER TIME UNIT
## per genome
head(get.mu(res))

## per nucleotide
mu <- get.mu(res, genome.size=1e4)
head(mu)

summary(mu)
hist(mu, border="lightgrey", col="grey", xlab="Mutation per day and nucleotide",
     main="Posterior distribution of mutation rate")

detach(fakeOutbreak)





cleanEx()
nameEx("simOutbreak")
### * simOutbreak

flush(stderr()); flush(stdout())

### Name: simple outbreak simulator
### Title: Simulation of pathogen genotypes during disease outbreaks
### Aliases: simOutbreak print.simOutbreak [.simOutbreak labels.simOutbreak
###   simOutbreak-class as.igraph.simOutbreak plot.simOutbreak disperse

### ** Examples

## Not run: 
##D dat <- list(n=0)
##D 
##D ## simulate data with at least 30 cases
##D while(dat$n < 30){
##D    dat <- simOutbreak(R0 = 2, infec.curve = c(0, 1, 1, 1), n.hosts = 100)
##D }
##D dat
##D 
##D ## plot first 30 cases
##D N <- dat$n
##D plot(dat[1:(min(N,30))], main="First 30 cases")
##D mtext(side=3, text="nb mutations / nb generations")
##D 
##D ## plot a random subset (n=10) of the first cases
##D x <- dat[sample(1:min(N,30), 10, replace=FALSE)]
##D plot(x, main="Random sample of 10 of the first 30 cases")
##D mtext(side=3, text="nb mutations / nb generations")
##D 
##D ## plot population dynamics
##D head(dat$dynam,15)
##D matplot(dat$dynam[1:max(dat$onset),],xlab="time",
##D    ylab="nb of individuals", pch=c("S","I","R"), type="b")
##D 
##D 
##D ## spatial model
##D w <-  exp(-sqrt((1:40)))
##D x <- simOutbreak(2, w, spatial=TRUE,
##D                  duration=500, disp=0.1, reach=.2)
##D 
##D ## spatial model, no dispersal
##D x <- simOutbreak(.5, w, spatial=TRUE,
##D                  duration=500, disp=0, reach=5)
## End(Not run)



cleanEx()
nameEx("tTree")
### * tTree

flush(stderr()); flush(stdout())

### Name: consensus ancestries
### Title: Simple transmission tree from outreaber's output
### Aliases: tTree get.tTree plot.tTree as.igraph.tTree findMutations.tTree
### Keywords: classes

### ** Examples

data(fakeOutbreak)
attach(fakeOutbreak)

## represent posterior ancestries
if(require(adegenet)){
transGraph(res, annot="", main="Posterior ancestries - support > 0.01",
   threshold=0.01, col.pal=spectral)
}
## get consensus ancestries
tre <- get.tTree(res)
plot(tre, annot="", main="Consensus ancestries")

## show match data/consensus ancestries
col <- rep("lightgrey", 30)
col[which(dat$ances != tre$ances)] <- "pink"
plot(tre, annot="", vertex.color=col, main="Consensus ancestries")
mtext(side=3, text="cases with erroneous ancestries in pink")


detach(fakeOutbreak)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
