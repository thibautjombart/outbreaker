
################################
## FUNCTION MAKING SIMULATIONS:
## simulates data, runs outbreaker, derive results, write files
##
makeSimul <- function(N=1, type=gsub(".*/","",getwd()){
    ## CHECKS ##
    if(!require(outbreaker)) stop("outbreaker is needed")
    if(!require(adegenet)) stop("adegenet is needed")
    if(!require(ape)) stop("ape is needed")
    if(!require(digest)) stop("digest is needed")
    if(!require(EpiEstim)) stop("EpiEstim is needed")



    ## AUXILIARY FUNCTIONS ##
      w.dens.gen <- function(mu, sigma, max.k){
        res <- sapply(0:max.k, DiscrSI, mu, sigma)
        names(res) <- 0:max.k
        return(res)
    }



    ## DEFINE PARAMETERS ##
    ## get type of change
    type.head <- unlist(strsplit(type,"-"))[1]

    ## R0
    if(type.head=="R"){
        R0 <- as.numeric(unlist(strsplit(type,"-"))[2])
    } else {
        R0 <- 1.5
    }

    ## generation time function (w)
    if(type.head=="w"){
        temp <- unlist(strsplit(type,"-"))[2]
        temp <- as.numeric(unlist(strsplit(temp, "_")))
        w.mu <- temp[1]
        w.sigma <- temp[2]
        w.k <- temp[3]
    } else {
        w.mu <- 2
        w.sigma <- 0.7
        w.k <- 10
    }

    w <- w.dens.gen(w.mu, w.sigma, w.k)

    ## mutation rates
    if(type.head=="mu1"){
       mu1 <- as.numeric(unlist(strsplit(type,"-"))[2])
       mu2 <- 0.25*mu1
    } else {
        mu1 <- 0.8e-4
        mu2 <- 0.25*mu1
    }

    ## imported cases
     if(type.head=="imp"){
       r.imp <- as.numeric(unlist(strsplit(type,"-"))[2])
    } else {
        r.imp <- 0.05
    }

    ## proportion of observed cases
    if(type.head=="samp"){
        p.samp <- as.numeric(unlist(strsplit(type,"-"))[2])
    } else {
        p.samp <- 1
    }

    ## proportion of observed cases
    if(type.head=="seq"){
        p.seq <- as.numeric(unlist(strsplit(type,"-"))[2])
    } else {
        p.seq <- 1
    }



    ## MAKE N SIMULATIONS ##
    for(i in 1:N){
        ## simulate outbreak ##
        full <- simOutbreak(R0=R0, infec.curve=w, mu.transi=mu1, mu.transv=mu2, rate.import.case=r.imp,
                            duration=100, n.hosts=200, seq.length=1e4, diverg.import=10)
        while(full$n < 10){
            full <- simOutbreak(R0=R0, infec.curve=w, mu.transi=mu1, mu.transv=mu2, rate.import.case=r.imp,
                            duration=100, n.hosts=200, seq.length=1e4, diverg.import=10)
        }

        ## subset data ##
        n.dat <- round(p.samp*full$n)
        dat <- full[sort(sample(1:n.dat))]

        ## check which cases are sequenced ##
        temp <- as.logical(rbinom(dat$n, prob=p.seq, size=1))
        which.seq <- which(temp)
        dna <- dat$dna[which.seq,,drop=FALSE]

        ## get collaction dates ##
        collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)

        ## run outbreaker ##
        BURNIN <- 2e4
        timing <- system.time(res <- outbreaker(dna=dna, dates=collecDates, idx.dna=which.seq, w.dens=w, init.tree="star", init.kappa=NULL,
                                      n.iter=1e5, sample.every=500,tune.every=500,burnin=BURNIN,find.import=TRUE, find.import.n=50,
                                      pi.param1=1, pi.param2=1, init.mu1=1e-5, init.gamma=1, move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE))

        ## extract information from results ##
        tre <- get.TTree.simple(res, burn=BURNIN)
        chains <- res$chains[res$chains$step>BURNIN,,drop=FALSE]

        stat <- list()
        stat$n <- dat$n

        ## prop of true ancestry in consensus
        stat$prop.ances.ok <- mean(tre$ances==dat$ances, na.rm=TRUE)

        ## mean support of actual ancestor
        alpha <- chains[,grep("alpha", names(chains))]
        alpha[is.na(alpha)] <- -1
        temp <- sapply(1:dat$n, function(i) mean(dat$ances[i]==alpha[,i],na.rm=TRUE))
        stat$msup.ances <- mean(temp, na.rm=TRUE)

        ## mean support of actual kappa
        kappa <- chains[,grep("kappa", names(chains))]
        temp <- sapply(1:dat$n, function(i) mean(dat$ngen[i]==kappa[,i],na.rm=TRUE))
        stat$msup.kappa <-mean(temp, na.rm=TRUE)

        ## mean error date
        Tinf <- chains[,grep("Tinf", names(chains))]
        temp <- abs(t(Tinf) - dat$dates)
        stat$merr.date <- mean(temp, na.rm=TRUE)

        ## proportion of successfully detected imported cases
        stat$prop.imp.ok <- sum(is.na(tre$ances[-1]) & is.na(dat$ances[-1]))/sum(is.na(dat$ances[-1]))

        ## mu1.ok: 1 if true value of mu1 is in the 95% CI, 0 otherwise
        stat$mu1.ok <- as.numeric(mu1>quantile(chains$mu1, 0.05) & mu1<quantile(chains$mu1, 0.95))

        ## mu2.ok: 1 if true value of mu2 is in the 95% CI, 0 otherwise
        stat$mu2.ok <- as.numeric(mu2>quantile(chains$mu2, 0.05) & mu2<quantile(chains$mu2, 0.95))

        ## merr.pi: the mean error in posterior values of pi
        stat$merr.pi <- mean(abs(chains$pi - p.samp))

        ## timing
        stat$time <- round(timing[3])

        ## hash key
        key <- paste(unlist(strsplit(digest(dat),""))[1:10], collapse="")
        stat <- data.frame(stat)
        row.names(stat) <- key

        ## SAVE OBJECTS TO FILE ##
        save(full,dat,collecDates,res,chains,tre,stat,key, file=paste(key,"Rdata",sep="."))



} # end makeSimul

## DATA SIMULATION ##

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)

## full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4, mu.transv=0.2e-4)
## dat <- full[1:20]
## collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
##save(w, full, dat, collecDates, file="Robjects/data4.RData")

##load("Robjects/data5.RData")
##plot(dat, main="Data")
############################################

BURNIN <- 2e4

w <- c(0,.1,.2,.5,2,.5,.2,.1)
full <- simOutbreak(R0=2, infec.curve=w, mu.transi=1e-4, mu.transv=0.2e-4)
dat <- full[1:50]
collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)
plot(dat, main="Data")



## ## TRY SPECIFYING THE RIGHT OUTLIERS, SEE IF IT WORKS ##
## initree <- rep(1, dat$n)
## initree[is.na(dat$ances)] <- 0
## moveances <- !is.na(dat$ances)

## system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree=initree, n.iter=2e5, move.ances=moveances))




############################################
## ESTIMATE EVERYTHING - PARALLEL VERSION ##
## run outbreaker
## system.time(res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=1e5))

system.time(res <- outbreaker.parallel(n.runs=4, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="star", n.iter=1e5))

## check results ##
plot.chains(res)

## check ancestries
x <- get.TTree.simple(res, burn=5e4)
temp <- x
temp$ances[is.na(temp$ances)] <- 0
temp2 <- dat$ances
temp2[is.na(temp2)] <- 0
temp2[1] <- NA
mean(temp$ances==temp2,na.rm=TRUE)

v.col <- rep("lightblue",length(x$ances))
notOk <- which(temp$ances!=temp2)
if(length(notOk)>0) v.col[notOk] <- "red"
par(mfrow=c(1,1))
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col)


## check frequency of external infections
alpha <- res$chains[,grep("alpha", names(res$chains))]
barplot(apply(alpha,2, function(e) mean(e==0)), main="Proportion of inferred external case")


## check mutation rates
par(mfrow=c(2,1))
plot.chains(res, "mu1",type="dens", omit=2e4)
abline(v=1e-4, col="blue")
plot.chains(res, "mu2",type="dens", omit=2e4)
abline(v=1e-4, col="blue")



## check kappa
v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$n.gen!=1)
if(length(notOk)>0) v.col[notOk] <- "red"
plot(dat,main="data", vertex.color=v.col)
x11();
plot(x,main="reconstruction (red=wrong kappa)", vertex.color=v.col, annot="n.gen")



## check Tinf
toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
boxplot(Tinf, col="grey")
points(dat$dates, col="red", pch="x",cex=2)



## check Pi
plot.chains(res, "pi", omit=BURNIN, type="de")
abline(v=1)



## check Phi
plot.chains(res, "phi", omit=BURNIN, type="de")
abline(v=1/20)

############################################
