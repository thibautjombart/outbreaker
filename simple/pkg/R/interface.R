
##################
## main functions
##################
outbreaker <- function(dna=NULL, dates, idx.dna=NULL,
                       w.dens, w.trunc=length(w.dens),
                       init.tree=c("seqTrack","random","star","none"),
                       init.kappa=NULL,
                       n.iter=1e5, sample.every=500, tune.every=500,
                       burnin=2e4, find.import=TRUE, find.import.n=50,
                       pi.param1=10, pi.param2=1,
                       init.mu1=NULL, init.gamma=1,
                       move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                       move.Tinf=TRUE, move.pi=TRUE,
                       quiet=TRUE, res.file.name="chains.txt", tune.file.name="tuning.txt", seed=NULL){

    ## CHECKS ##
    if(!require(ape)) stop("the ape package is required but not installed")

    ## HANDLE MISSING DNA ##
    useDna <- !is.null(dna)
    if(is.null(dna)){
        dna <- as.DNAbin(matrix('a',ncol=1,nrow=length(dates)))
        move.mut <- FALSE
        find.import <- FALSE
        init.tree <- "star"
    }

    if(!inherits(dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if(!is.matrix(dna)) dna <- as.matrix(dna)
    if(is.character(dates)) stop("dates are characters; they must be integers or dates with POSIXct format (see ?as.POSIXct)")
    if(is.character(init.tree)) {
        init.tree <- match.arg(init.tree)
    } else {
        if(length(init.tree) != length(dates)) stop("inconvenient length for init.tree")
        init.tree[is.na(init.tree)|init.tree<1] <- 0
        if(max(init.tree)>length(dates)) stop("inconvenient values in init.tree (some indices > n)")
        ances <- as.integer(init.tree-1) # translate indices on C scale (0:(n-1))
    }
    if(length(w.dens)<w.trunc) stop(paste("incomplete w.dens: values needed from t=0 to t=", w.trunc-1,sep=""))
    w.dens[1] <- 0 # force w_0 = 0
    w.dens[w.dens<0] <- 0
    if(sum(w.dens) < 1e-14) stop("w.dens is zero everywhere")
    if(init.mu1<0) stop("init.mu1 < 0")
    if(init.gamma<0) stop("init.gamma < 0")



    ## PROCESS INPUTS ##
    ## dna ##
    n.seq <- as.integer(nrow(dna))
    n.ind <- as.integer(length(dates))
    n.nucl <- as.integer(ncol(dna))
    dnaraw <- unlist(as.list(dna),use.names=FALSE)
    ## if(n.ind != length(dates)) stop(paste("dna and dates have different number of individuals -",n.ind,"versus",length(dates)))

    ## handle dates ##
    if(is.numeric(dates)){
        if(sum(abs(dates-round(dates))>1e-15)) warning("dates have been rounded to nearest integers")
        dates <- as.integer(round(dates))
    }

    if(inherits(dates, "POSIXct")){
        dates <- difftime(dates, min(dates), units="days")
    }
    dates <- as.integer(dates)


    ## handle idx.dna ##
    ## need to go from: id of case for each sequence
    ## to: position of the sequence in DNA matrix for each case
    ## -1 is used for missing sequences
    if(is.null(idx.dna)) {
        idx.dna <- 1:nrow(dna)
    }

    if(any(!idx.dna %in% 1:n.ind)) stop("DNA sequences provided for unknown cases (some idx.dna not in 1:n.ind)")
    if(length(idx.dna)!=nrow(dna)) stop("length of idx.dna does not match the number of DNA sequences")

    idx.dna.for.cases <- match(1:n.ind, idx.dna)
    idx.dna.for.cases[is.na(idx.dna.for.cases)] <- 0
    idx.dna.for.cases <- as.integer(idx.dna.for.cases-1) # for C

    ## check generation time function ##
    w.dens <- as.double(w.dens)
    w.dens <- w.dens/sum(w.dens)
    if(any(is.na(w.dens))) stop("NAs in w.dens after normalization")
    w.trunc <- as.integer(w.trunc)

    ## init.kappa ##
    ## if NULL, will be ML assigned (code is kappa_i<0)
    if(is.null(init.kappa)) init.kappa <- rep(0L,n.ind)
    init.kappa <- as.integer(rep(init.kappa, length=n.ind)) #recycle


    ## find initial tree ##
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    canBeAnces <- outer(dates,dates,FUN="<") # strict < is needed as we impose w(0)=0
    diag(canBeAnces) <- FALSE

    if(is.character(init.tree)){
        ## seqTrack init
        if(init.tree=="seqTrack"){
            if(!all(1:n.ind %in% idx.dna)) stop("can't use seqTrack initialization with missing DNA sequences")
            D <- as.matrix(dist.dna(dna, model="TN93"))
            D[!canBeAnces] <- 1e15
            ances <- apply(D,2,which.min)-1 # -1 for compatibility with C
            ances[dates==min(dates)] <- -1 # unknown ancestor
            ances <- as.integer(ances)
        }

        ## random init
        if(init.tree=="random"){
            ances <- apply(canBeAnces, 2, function(e) ifelse(length(which(e))>0, sample(which(e),1), NA) )
            ances <- ances-1
            ances[is.na(ances)] <- -1L
            ances <- as.integer(ances)
        }

        ## star-shape init
        if(init.tree=="star"){
            ances <- rep(which.min(dates), length(dates))
            ances[dates==min(dates)] <- 0
            ances <- as.integer(ances-1) # put on C scale
        }

        ## no ancestry init
        if(init.tree=="none"){
            ances <- as.integer(rep(-1,length(dates)))
        }
    }

    ## handle seed ##
    if(is.null(seed)){
        seed <- as.integer(runif(1,min=0,max=2e9))
    }

    ## handle find.import ##
    if(find.import){
        find.import.n <- max(find.import.n,30) # import at least based on 30 values
        find.import.at <- as.integer(round(burnin + find.import.n*sample.every))
        if(find.import.at>=n.iter) stop(paste("n.iter (", n.iter, ") is less than find.import.at (", find.import.at,")", sep=""))
    } else {
        find.import.at <- as.integer(0)
    }

    ## coerce type for remaining arguments ##
    n.iter <- as.integer(n.iter)
    sample.every <- as.integer(sample.every)
    tune.every <- as.integer(tune.every)
    pi.param1 <- as.double(pi.param1)
    pi.param2 <- as.double(pi.param2)
    phi.param1 <- as.double(1)
    phi.param2 <- as.double(10)
    if(is.null(init.mu1)) {
        init.mu1 <- 0.5/ncol(dna)
    }
    init.mu1 <- as.double(init.mu1)
    init.gamma <- as.double(init.gamma)
    move.mut <- as.integer(move.mut)
    move.ances <- as.integer(rep(move.ances, length=n.ind))
    move.kappa <- as.integer(rep(move.kappa, length=n.ind))
    move.Tinf <- as.integer(move.Tinf)
    move.pi <- as.integer(move.pi)
    quiet <- as.integer(quiet)
    res.file.name <- as.character(res.file.name)[1]
    tune.file.name <- as.character(tune.file.name)[1]
    burnin <- as.integer(burnin)
    find.import.int <- as.integer(find.import)

    ## create empty output vector for genetic distances ##
    dna.dist <- integer(n.ind*(n.ind-1)/2)
    stopTuneAt <- integer(1)

    temp <- .C("R_outbreaker",
               dnaraw, dates, n.ind, n.seq, n.nucl,  idx.dna.for.cases,
               w.dens, w.trunc,
               ances, init.kappa, n.iter, sample.every, tune.every,
               pi.param1, pi.param2, init.mu1, init.gamma,
               move.mut, move.ances, move.kappa, move.Tinf,
               move.pi,
               find.import.int, burnin, find.import.at, quiet,
               dna.dist, stopTuneAt, res.file.name, tune.file.name, seed,
               PACKAGE="outbreaker")

    D <- temp[[27]]
    D[D<0] <- NA
    stopTuneAt <- temp[[28]]

    cat("\nComputations finished.\n\n")

    ## make D a 'dist' object ##
    attr(D,"Size") <- n.ind
    attr(D,"Diag") <- FALSE
    attr(D,"Upper") <- FALSE
    class(D) <- "dist"


    ## BUILD OUTPUT ##
    chains <- read.table(res.file.name, header=TRUE)
    chains$run <- rep(1, nrow(chains))
    call <- match.call()
    res <- list(chains=chains, collec.dates=dates, w=w.dens[1:w.trunc], D=D, idx.dna=idx.dna, tune.end=stopTuneAt,
                find.import=find.import, burnin=burnin, find.import.at=find.import.at, n.runs=1, call=call)

    return(res)
} # end outbreaker









###############################
## version with multiple runs
###############################
outbreaker.parallel <- function(n.runs, parallel=require("parallel"), n.cores=NULL,
                                dna=NULL, dates, idx.dna=NULL, w.dens, w.trunc=length(w.dens),
                                init.tree=c("seqTrack","random","star","none"),
                                init.kappa=NULL,
                                n.iter=1e5, sample.every=500, tune.every=500,
                                burnin=2e4, find.import=TRUE, find.import.n=50,
                                pi.param1=10, pi.param2=1,
                                init.mu1=NULL, init.gamma=1,
                                move.mut=TRUE, move.ances=TRUE, move.kappa=TRUE,
                                move.Tinf=TRUE, move.pi=TRUE, move.phi=TRUE,
                                quiet=TRUE, res.file.name="chains.txt", tune.file.name="tuning.txt", seed=NULL){

    ## SOME CHECKS ##
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel:::detectCores()
    }


    ## GET FILE NAMES ##
    res.file.names <- paste("run", 1:n.runs, "-", res.file.name, sep="")
    tune.file.names <- paste("run", 1:n.runs, "-", tune.file.name, sep="")


    ## HANDLE SEED ##
    if(is.null(seed)){
        seed <- as.integer(runif(n.runs,min=0,max=2e9))
    } else {
        seed <- rep(seed, length=n.runs)
    }


    ## COMPUTATIONS ##
    if(parallel){
        ## create cluster ##
        clust <- makeCluster(n.cores)

        ## load outbreaker for each child ##
        clusterEvalQ(clust, library(outbreaker))

        ## transfer data onto each child ##
        listArgs <- c("dna", "dates", "idx.dna", "w.dens", "w.trunc", "init.tree", "init.kappa", "n.iter", "sample.every", "tune.every", "burnin", "find.import", "find.import.n", "pi.param1", "pi.param2", "init.mu1", "init.gamma", "move.mut", "move.ances", "move.kappa", "move.Tinf", "move.pi", "move.phi", "res.file.names", "tune.file.names", "seed")

        clusterExport(clust, listArgs, envir=environment())

        ## set calls to outbreaker on each child ##
        res <- parLapply(clust, 1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna, w.dens=w.dens, w.trunc=w.trunc,
                                                           init.tree=init.tree, init.kappa=init.kappa,
                                                           n.iter=n.iter, sample.every=sample.every,
                                                           tune.every=tune.every, burnin=burnin,
                                                           find.import=find.import, find.import.n=find.import.n,
                                                           pi.param1=pi.param1, pi.param2=pi.param2,
                                                           init.mu1=init.mu1, init.gamma=init.gamma,
                                                           move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
                                                           move.Tinf=move.Tinf, move.pi=move.pi,
                                                           quiet=TRUE, res.file.name=res.file.names[i],
                                                           tune.file.name=tune.file.names[i], seed=seed[i]))

        ## close parallel processes ##
        stopCluster(clust)

        ## Version with mclapply - doesn't work on windows ##
        ## res <- mclapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna, w.dens=w.dens, w.trunc=w.trunc,
        ##                                                     init.tree=init.tree, init.kappa=init.kappa,
        ##                                                     n.iter=n.iter, sample.every=sample.every,
        ##                                                     tune.every=tune.every, burnin=burnin,
        ##                                                     find.import=find.import, find.import.n=find.import.n,
        ##                                                     pi.param1=pi.param1, pi.param2=pi.param2,
        ##                                                     init.mu1=init.mu1, init.gamma=init.gamma,
        ##                                                     move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
        ##                                                     move.Tinf=move.Tinf, move.pi=move.pi,
        ##                                                     quiet=TRUE, res.file.name=res.file.names[i],
        ##                                                     tune.file.name=tune.file.names[i], seed=seed[i]),
        ##                   mc.cores=n.cores, mc.silent=FALSE, mc.cleanup=TRUE, mc.preschedule=TRUE, mc.set.seed=TRUE)
    } else {
        res <- lapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, idx.dna=idx.dna, w.dens=w.dens, w.trunc=w.trunc,
                                                        init.tree=init.tree, init.kappa=init.kappa,
                                                        n.iter=n.iter, sample.every=sample.every,
                                                        tune.every=tune.every, burnin=burnin,
                                                        find.import=find.import, find.import.n=find.import.n,
                                                        pi.param1=pi.param1, pi.param2=pi.param2,
                                                        init.mu1=init.mu1, init.gamma=init.gamma,
                                                        move.mut=move.mut, move.ances=move.ances, move.kappa=move.kappa,
                                                        move.Tinf=move.Tinf, move.pi=move.pi,
                                                        quiet=TRUE, res.file.name=res.file.names[i],
                                                        tune.file.name=tune.file.names[i], seed=seed[i]))
    }


    ## MERGE RESULTS ##
    res.old <- res
    res <- res[[1]]
    res$tune.end <- max(sapply(res.old, function(e) e$tune.end))
    res$chains <- Reduce(rbind, lapply(res.old, function(e) e$chains))
    res$chains$run <- rep(1:n.runs, each=nrow(res.old[[1]]$chains))
    res$n.runs <- n.runs
    res$call <- match.call()

    ## RETURN RESULTS ##
    return(res)
} # end outbreaker.parallel

