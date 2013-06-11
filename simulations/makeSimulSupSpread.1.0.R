
################################
## FUNCTION MAKING SIMULATIONS:
## simulates data, runs outbreaker, derive results, write files
##
makeSimulZuper <- function(N=1){
    ## CHECKS ##
    if(!require(outbreaker)) stop("outbreaker is needed")
    if(!require(adegenet)) stop("adegenet is needed")
    if(!require(ape)) stop("ape is needed")
    if(!require(digest)) stop("digest is needed")
    if(!require(EpiEstim)) stop("EpiEstim is needed")

    ## create output directory if needed
    if(length(dir(pattern="superSpread"))==0) dir.create("superSpread")
    odir <- getwd()
    setwd("superSpread")
    on.exit(setwd(odir))

    ## AUXILIARY FUNCTIONS ##
    w.dens.gen <- function(mu, sigma, max.k){
        res <- sapply(0:max.k, DiscrSI, mu, sigma)
        names(res) <- 0:max.k
        return(res)
    }


    ## DEFINE PARAMETERS ##
    ## generation time distribution (from Cori et al., submitted)
    w.mu <-  8
    w.sigma <- 4
    w.k <- 100
    w <- w.dens.gen(w.mu, w.sigma, w.k)

    ## mutation rates
    mu1 <- 1.5e-4
    mu2 <- mu1/4

    ## basic reproduction numbers
    R0 <- c(1.5,20)

    ## frequency of super-spreaders
    f.grp <- c(0.95,0.05)


    ## MAKE N SIMULATIONS ##
    for(i in 1:N){
        ## simulate outbreak ##
        dat <- simOutbreak(R0=R0, infec.curve=w, mu.transi=mu1, mu.transv=mu2, rate.import.case=0,
                            duration=100, n.hosts=100, seq.length=1e4, diverg.import=10, group.freq=f.grp)

        ## run simulations until at least 10 cases and 1 super spreader
        noZuper <- TRUE
        while(dat$n < 10 && sum(dat$group==2)<1 && noZuper){
            dat <- simOutbreak(R0=R0, infec.curve=w, mu.transi=mu1, mu.transv=mu2, rate.import.case=0.0,
                            duration=100, n.hosts=100, seq.length=1e4, diverg.import=10, group.freq=f.grp)

            ## check if there is at least one super-spreader
            temp <- sapply(1:dat$n, function(i) sum(dat$ances %in% i))
            noZuper <- !any(dat$group==2 & temp>=quantile(temp,.95))
        }


        ## get collection dates ##
        collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)

        ## run outbreaker ##
        BURNIN <- 2e4

        ## hash key
        key <- paste(unlist(strsplit(digest(dat),""))[1:10], collapse="")

        ## create dir and move to it
        dir.create(key)
        curdir <- getwd()
        setwd(key)

        cat("Ongoing computations", file="ONGOING")

        ## run outbreaker - know that all outbreak sampled
        res <- outbreaker(dna=dat$dna, dates=collecDates, w.dens=w, init.tree="seqTrack", init.kappa=1,
                          n.iter=1e5, sample.every=500,tune.every=500,burnin=BURNIN,find.import=FALSE,
                          pi.param1=1, pi.param2=1, init.mu1=1e-5, init.gamma=1,
                          move.mut=TRUE, move.ances=TRUE, move.kappa=FALSE)

        res.nodna <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, init.tree="seqTrack", init.kappa=1,
                           n.iter=1e5, sample.every=500,tune.every=500,burnin=BURNIN,find.import=FALSE,
                           pi.param1=1, pi.param2=1, init.mu1=1e-5, init.gamma=1,
                           move.mut=TRUE, move.ances=TRUE, move.kappa=FALSE)

        ## extract information from results ##
        tre <- get.TTree.simple(res, burn=BURNIN)
        tre.nodna <- get.TTree.simple(res.nodna, burn=BURNIN)
        chains <- res$chains[res$chains$step>BURNIN,,drop=FALSE]

        stat <- list()
        stat$type <- "zuper"
        stat$n <- dat$n

        ## prop of true ancestry in consensus
        stat$prop.ances.ok <- mean(tre$ances==dat$ances, na.rm=TRUE)

        ## mean support of actual ancestor
        alpha <- chains[,grep("alpha", names(chains))]
        alpha[is.na(alpha)] <- -1
        temp <- sapply(1:dat$n, function(i) mean(dat$ances[i]==alpha[,i],na.rm=TRUE))
        stat$msup.ances <- mean(temp, na.rm=TRUE)

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
        stat$merr.pi <- abs(mean(chains$pi) - 1)

        ## turn stat into a data.frame
        stat <- data.frame(stat)
        row.names(stat) <- key


        ## SUPER-SPREADER SPECIFIC ANALYSES ##
        ## refine groups - only realized super-spreader are 'zuper'
        temp <- sapply(1:dat$n, function(i) sum(dat$ances %in% i))
        trueZuper <- which(dat$group==2 & temp>=quantile(temp,.95))

        ## check if there is at least one super-spreader
        if(length(trueZuper)>0){
            group <- rep(1,dat$n)
            group[trueZuper] <- 2
            group.adjfreq <- as.numeric(table(group)/dat$n)

            ## chi2 test
            table(group)
            table(group[tre$ances])
            chitest <- chisq.test(table(group[tre$ances]), p=group.adjfreq, simulate=TRUE)
            chitest$p.value
            chitest.nodna <- chisq.test(table(group[tre.nodna$ances]), p=group.adjfreq, simulate=TRUE)
            chitest.nodna$p.value

            ## distribution of number of descendents
            nbDesc <- sapply(1:dat$n, function(i) sum(tre$ances %in% i))
            names(nbDesc) <-  1:dat$n
            nbDesc.nodna <- sapply(1:dat$n, function(i) sum(tre.nodna$ances %in% i))
            names(nbDesc.nodna) <-  1:dat$n

            ## see which one is identified as superspreader
            areFoundZuper <- nbDesc[group==2] >= quantile(nbDesc, 0.95)
            areFoundZuper.nodna <- nbDesc.nodna[group==2] >= quantile(nbDesc.nodna, 0.95)

            ## output to file
            cat(c(mean(areFoundZuper),mean(areFoundZuper.nodna)), file=paste(key,".zuper.txt", sep=""))
        } else {
            cat(NA, file=paste(key,".zuper.txt", sep=""))
            areFoundZuper <- NULL
            areFoundZuper.nodna <- NULL
            nbDesc <- NULL
            nbDesc.nodna <- NULL
            group <- NULL
            chitest <- NULL
        }


        ## SAVE OBJECTS TO FILE ##
        save(dat,collecDates,res,res.nodna,chains,tre,tre.nodna,stat,group, chitest, nbDesc, key, areFoundZuper, areFoundZuper.nodna, file=paste(key,"RData",sep="."))

        ## ## COMPILE/WRITE INPUT DATA ##
        ## inputs <- list()
        ## inputs$R <- R0
        ## inputs$w.mu <- w.mu
        ## inputs$w.sigma <- w.sigma
        ## inputs$w.k <- w.k
        ## inputs$mu1 <- mu1
        ## inputs$mu2 <- mu2
        ## inputs$r.imp.cases <- 0
        ## inputs$p.obs.cases <- 1
        ## inputs$p.seq.data <- 1

        ## inputs <- data.frame(inputs)
        ## row.names(inputs) <- key

        ## write.csv(inputs, file=paste(key,".in.csv", sep=""))


        ## WRITE STAT TO FILE ##
        write.csv(stat, file=paste(key,".out.csv", sep=""))

        ## GO BACK TO PREVIOUS WORKING DIRECTORY ##
        file.remove("ONGOING")
        setwd(curdir)

    } # end for loop

} # end makeSimulZuper
