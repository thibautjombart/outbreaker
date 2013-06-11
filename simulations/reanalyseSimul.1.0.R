
################################
## FUNCTION MAKING SIMULATIONS:
## simulates data, runs outbreaker, derive results, write files
##
reanalyseSimul <- function(key){
    ## CHECKS ##
    if(!require(outbreaker)) stop("outbreaker is needed")
    if(!require(adegenet)) stop("adegenet is needed")
    if(!require(ape)) stop("ape is needed")
    if(!require(digest)) stop("digest is needed")
    if(!require(EpiEstim)) stop("EpiEstim is needed")

    ## EXIT IF ONGOING ANALYSES ##
    path <- key
    ## if(length(dir(path, pattern="ONGOING"))>0){
    ##     cat(paste("\nAnalyses ongoing in directory",path,"\n"))
    ##     return()
    ## }


    ## REDO ANALYSIS WITH FIXED MUTATION RATES ##
    ## GET INTO DIRECTORY
    odir <- getwd()
    setwd(path)
    on.exit(setwd(odir))
    cat("Ongoing computations", file="ONGOING")

    ## LOAD DATA
    BURNIN <- 2e4
    ##    load(dir(pattern="RData"))
    load(paste(key, ".RData", sep=""))

    ## GET PROPORTION OF SEQUENCED CASES
    p.seq <- read.csv(paste(key,".in.csv",sep=""))$p.seq.data
    temp <- as.logical(rbinom(dat$n, prob=p.seq, size=1))
    which.seq <- which(temp)
    dna <- dat$dna[which.seq,,drop=FALSE]

    ## GET MUTATION RATES
    mu1 <- read.csv(paste(key,".in.csv",sep=""))$mu1
    if(mu1<2e-12){
        gamma <- 1.0
    } else {
        gamma <- read.csv(paste(key,".in.csv",sep=""))$mu2/mu1
    }

    ## run outbreaker ##
    BURNIN <- 2e4
    timing <- system.time(res2 <- outbreaker(dna=dna, dates=collecDates, idx.dna=which.seq, w.dens=res$w, init.tree="star", init.kappa=NULL,
                                                n.iter=1e5, sample.every=500,tune.every=500,burnin=BURNIN,find.import=TRUE, find.import.n=50,
                                                pi.param1=1, pi.param2=1, init.mu1=mu1, init.gamma=gamma, move.mut=FALSE, move.ances=TRUE, move.kappa=TRUE))

    ## extract information from results ##
    tre2 <- get.TTree.simple(res2, burn=BURNIN)
    chains <- res2$chains[res2$chains$step>BURNIN,,drop=FALSE]

    stat2 <- list()
    stat2$type <- paste(stat$type,"fixedMu",sep="-")
    stat2$n <- dat$n

    ## prop of true ancestry in consensus
    stat2$prop.ances.ok <- mean(tre2$ances==dat$ances, na.rm=TRUE)

    ## mean support of actual ancestor
    alpha <- chains[,grep("alpha", names(chains))]
    alpha[is.na(alpha)] <- -1
    temp <- sapply(1:dat$n, function(i) mean(dat$ances[i]==alpha[,i],na.rm=TRUE))
    stat2$msup.ances <- mean(temp, na.rm=TRUE)

    ## mean support of actual kappa
    kappa <- chains[,grep("kappa", names(chains))]
    temp <- sapply(1:dat$n, function(i) mean(dat$ngen[i]==kappa[,i],na.rm=TRUE))
    stat2$msup.kappa <-mean(temp, na.rm=TRUE)

    ## mean error date
    Tinf <- chains[,grep("Tinf", names(chains))]
    temp <- abs(t(Tinf) - dat$dates)
    stat2$merr.date <- mean(temp, na.rm=TRUE)

    ## proportion of successfully detected imported cases
    stat2$prop.imp.ok <- sum(is.na(tre2$ances[-1]) & is.na(dat$ances[-1]))/sum(is.na(dat$ances[-1]))

    ## mu.ok: 1 if true value of mu (=mu1+mu2) is in the 95% CI, 0 otherwise
    stat2$mu.ok <- 1

    ## mu1.ok: 1 if true value of mu1 is in the 95% CI, 0 otherwise
    stat2$mu1.ok <- 1

    ## mu2.ok: 1 if true value of mu2 is in the 95% CI, 0 otherwise
    stat2$mu2.ok <- 1

    ## merr.pi: the mean error in posterior values of pi
    p.samp <- read.csv(paste(key,".in.csv",sep=""))$p.obs.cases
    stat2$merr.pi <- abs(mean(chains$pi) - p.samp)

    ## timing
    stat2$time <- round(timing[3])

    ## turn stat2 into a data.frame
    stat2 <- data.frame(stat2)
    row.names(stat2) <- key

    ## SAVE OBJECTS TO FILE ##
    save(res2,tre2,stat2,key, file=paste(key,"-2.RData",sep=""))


    ## COMPILE/WRITE INPUT DATA ##
    file.copy(paste(key,".in.csv",sep=""),paste(key,"-2.in.csv",sep=""))

    ## WRITE STAT TO FILE ##
    write.csv(stat2, file=paste(key,"-2.out.csv", sep=""))

    ## GO BACK TO PREVIOUS WORKING DIRECTORY ##
    file.remove("ONGOING")
    setwd(odir)

} # end makeSimul
