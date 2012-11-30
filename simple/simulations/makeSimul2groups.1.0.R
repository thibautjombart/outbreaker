
################################
## FUNCTION MAKING SIMULATIONS:
## simulates data, runs outbreaker, derive results, write files
##
makeSimul2groups <- function(N=1){
    ## CHECKS ##
    if(!require(outbreaker)) stop("outbreaker is needed")
    if(!require(adegenet)) stop("adegenet is needed")
    if(!require(ape)) stop("ape is needed")
    if(!require(digest)) stop("digest is needed")
    if(!require(EpiEstim)) stop("EpiEstim is needed")

    ## create output directory if needed
    if(length(dir(pattern="sim2groups"))==0) dir.create("sim2groups")
    odir <- getwd()
    setwd("sim2groups")
    on.exit(setwd(odir))

    ## AUXILIARY FUNCTIONS ##
    w.dens.gen <- function(mu, sigma, max.k){
        res <- sapply(0:max.k, DiscrSI, mu, sigma)
        names(res) <- 0:max.k
        return(res)
    }


    ## DEFINE PARAMETERS ##
    ## generation time distribution (from Cori et al., submitted)
    ## flu-like distribution
    w.mu <-  2.5
    w.sigma <- 1.5
    w.k <- 100
    w <- w.dens.gen(w.mu, w.sigma, w.k)

    ## mutation rates
    mu1 <- 1e-4
    mu2 <- mu1/4

    ## basic reproduction numbers
    R0 <- c(1.5,3)

    ## frequency of super-spreaders
    f.grp <- c(0.5,0.5)


    ## MAKE N SIMULATIONS ##
    for(i in 1:N){
        ## simulate outbreak ##
        dat <- simOutbreak(R0=R0, infec.curve=w, mu.transi=mu1, mu.transv=mu2, rate.import.case=0,
                           duration=100, n.hosts=100, seq.length=1e4, diverg.import=10, group.freq=f.grp)

        ## run simulations until at least 10 cases and significant differences between groups
        noDiff <- TRUE

        while(dat$n < 20 || noDiff){
            dat <- simOutbreak(R0=R0, infec.curve=w, mu.transi=mu1, mu.transv=mu2, rate.import.case=0,
                               duration=100, n.hosts=100, seq.length=1e4, diverg.import=10, group.freq=f.grp)

            ## check that both groups are there
            bothGrp <- all(table(dat$group)>5)

            ## check differences between groups
            if(bothGrp){
                noDiff <- chisq.test(table(dat$group[dat$ances]), p=c(.5,.5))$p.value>0.01
            } else {
                noDiff <- TRUE
            }

        }


        ## get collection dates ##
        collecDates <- dat$dates+sample(0:(length(w)-1), length(dat$dates), replace=TRUE, prob=w)

        ## run outbreaker ##
        BURNIN <- 2e4

        ## hash key
        key <- paste(unlist(strsplit(digest(dat),""))[1:10], collapse="")
        cat("\nsimulation: ", key)

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

        ## run outbreaker - know that all outbreak sampled
        res.nodna <- outbreaker(dna=NULL, dates=collecDates, w.dens=w, init.tree="seqTrack", init.kappa=1,
                                n.iter=1e5, sample.every=500,tune.every=500,burnin=BURNIN,find.import=FALSE,
                                pi.param1=1, pi.param2=1, init.mu1=1e-5, init.gamma=1,
                                move.mut=FALSE, move.ances=TRUE, move.kappa=FALSE)

        ## extract information from results ##
        tre <- get.TTree.simple(res, burn=BURNIN)
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


        ## GROUP-DIFFERENCES SPECIFIC ANALYSES ##
        ## test if there are actually differences
        temp <- sapply(1:dat$n, function(i) sum(i==dat$ances,na.rm=TRUE))

        ## check if there is at least one super-spreader
        if(t.test(temp~dat$group)$p.value>.01){
            ## get R distributions
            Rval <- get.R(res, burnin=BURNIN)
            Rval1 <- apply(Rval[,dat$group==1],1,mean)
            Rval2 <- apply(Rval[,dat$group==2],1,mean)

            Rval.nodna <- get.R(res.nodna, burnin=BURNIN)
            Rval1.nodna <- apply(Rval.nodna[,dat$group==1],1,mean)
            Rval2.nodna <- apply(Rval.nodna[,dat$group==2],1,mean)

            ## make tests
            myTest <- t.test(Rval1, Rval2)
            myTest.nodna <- t.test(Rval1.nodna, Rval2.nodna)

            ## output to file
            output <- data.frame(meanR1=mean(Rval1),meanR2=mean(Rval2), pval=myTest$p.value,
                                 meanR1.nodna=mean(Rval1.nodna),meanR2.nodna=mean(Rval2.nodna),
                                 pval.nodna=myTest.nodna$p.value)
            write.csv(output, file=paste(key,".test2groups.csv", sep=""))

            ## SAVE OBJECTS TO FILE ##
            save(dat, collecDates, res, res.nodna, chains, tre, stat, key, myTest, myTest.nodna, file=paste(key,"RData",sep="."))

            ## WRITE STAT TO FILE ##
            write.csv(stat, file=paste(key,".out.csv", sep=""))

            ## GO BACK TO PREVIOUS WORKING DIRECTORY ##
            file.remove("ONGOING")
            setwd(curdir)
            cat(paste("\n",key,"OK\n"))

        } else {
            ## if there was no significant difference in the data (not in the inference), remove and exit
            file.remove("ONGOING")
            cat(paste("\n",key,"removed\n"))

            setwd(curdir)
            cmd <- paste("rm",key,"-r &")
            system(cmd)
        }

    } # end for loop

} # end makeSimulZuper
