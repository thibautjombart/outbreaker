
################################
## FUNCTION MAKING SIMULATIONS:
## simulates data, runs outbreaker, derive results, write files
##
makeSimul <- function(N=1, type=gsub(".*/","",getwd())){
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

        ## hash key
        key <- paste(unlist(strsplit(digest(dat),""))[1:10], collapse="")

        ## create dir and move to it
        dir.create(key)
        odir <- getwd()
        setwd(key)

        cat("Ongoing computations", file="ONGOING")

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
        stat$merr.pi <- abs(mean(chains$pi) - p.samp)

        ## timing
        stat$time <- round(timing[3])

        ## turn stat into a data.frame
        stat <- data.frame(stat)
        row.names(stat) <- key

        ## SAVE OBJECTS TO FILE ##
        save(full,dat,collecDates,res,chains,tre,stat,key, file=paste(key,"Rdata",sep="."))


        ## COMPILE/WRITE INPUT DATA ##
        inputs <- list()
        inputs$R <- R0
        inputs$w.mu <- w.mu
        inputs$w.sigma <- w.sigma
        inputs$w.k <- w.k
        inputs$mu1 <- mu1
        inputs$mu2 <- mu2
        inputs$r.imp.cases<- r.imp
        inputs$p.obs.cases <- p.samp
        inputs$p.seq.data <- p.seq

        inputs <- data.frame(inputs)
        row.names(inputs) <- key

        write.csv(inputs, file=paste(key,".in.csv", sep=""))


        ## WRITE STAT TO FILE ##
        write.csv(stat, file=paste(key,".out.csv", sep=""))

        ## GO BACK TO PREVIOUS WORKING DIRECTORY ##
        file.remove("ONGOING")
        setwd(odir)

    } # end for loop

} # end makeSimul
