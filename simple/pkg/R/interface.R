
##################
## main functions
##################
outbreaker <- function(dna, dates, w.dens, w.trunc=length(w.dens),
                       init.tree=c("seqTrack","random","star","none"),
                       n.iter=2e6, sample.every=1000, tune.every=1000,
                       pi.param1=10, pi.param2=1, quiet=TRUE){
    ## CHECKS ##
    if(!require(ape)) stop("the ape package is required but not installed")
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


    ## PROCESS INPUTS ##
    ## dna ##
    n.ind <- as.integer(nrow(dna))
    n.nucl <- as.integer(ncol(dna))
    dnaraw <- unlist(as.list(dna),use.names=FALSE)
    if(n.ind != length(dates)) stop(paste("dna and dates have different number of individuals -",n.ind,"versus",length(dates)))

    ## dates in numeric ##
    if(is.numeric(dates)){
        if(sum(abs(dates-round(dates))>1e-15)) warning("dates have been rounded to nearest integers")
        dates <- as.integer(round(dates))
    }

    ## if dates are in POSIXct format ##
    if(inherits(dates, "POSIXct")){
        dates <- as.integer(difftime(dates, min(dates), units="days"))
    }

    ## make sure minimugit pum date is 0 ##
    dates <- as.integer(dates-min(dates))

    ## check generation time function ##
    w.dens <- as.double(w.dens)
    w.dens <- w.dens/sum(w.dens)
    if(any(is.na(w.dens))) stop("NAs in w.dens after normalization")
    w.trunc <- as.integer(w.trunc)


    ## find initial tree ##

    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    canBeAnces <- outer(dates,dates,FUN="<=")
    diag(canBeAnces) <- FALSE

    if(is.character(init.tree)){
        ## seqTrack init
        if(init.tree=="seqTrack"){
            D <- as.matrix(dist.dna(dna, model="TN93"))
            D[!canBeAnces] <- 1e15
            ances <- apply(D,2,which.min)-1 # -1 for compatibility with C
            ances[which.min(dates)] <- -1 # unknown ancestor
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

    ## coerce type for remaining arguments ##
    n.iter <- as.integer(n.iter)
    sample.every <- as.integer(sample.every)
    tune.every <- as.integer(tune.every)
    pi.param1 <- as.double(pi.param1)
    pi.param2 <- as.double(pi.param2)
    quiet <- as.integer(quiet)

    ## create empty output vector for genetic distances ##
    dna.dist <- integer(n.ind*(n.ind-1)/2)
    stopTuneAt <- integer(1)

    temp <- .C("R_outbreaker",
       dnaraw, dates, n.ind, n.nucl,
       w.dens, w.trunc,
       ances, n.iter, sample.every, tune.every, pi.param1, pi.param2, quiet,
       dna.dist, stopTuneAt,
       PACKAGE="outbreaker")

    D <- temp[[14]]
    stopTuneAt <- temp[[15]]

    cat("\nComputations finished.\n\n")

    ## make D a 'dist' object ##
    attr(D,"Size") <- n.ind
    attr(D,"Diag") <- FALSE
    attr(D,"Upper") <- FALSE
    class(D) <- "dist"


    ## BUILD OUTPUT ##
    chains <- read.table("output.txt",header=TRUE)
    call <- match.call()
    res <- list(chains=chains, collec.dates=dates, w=w.dens[1:w.trunc], D=D, tune.end=stopTuneAt, call=call)

    return(res)
} # end outbreaker









###############################
## version with multiple runs
###############################
outbreaker.parallel <- function(n.runs, multicore=require("multicore"), n.cores=NULL,
                                dna, dates, w.dens, w.trunc=length(w.dens),
                                init.tree=c("seqTrack","random","star"),
                                n.iter=2e6, sample.every=1000, tune.every=1000,
                                pi.param1=10, pi.param2=1, quiet=TRUE){

    ## SOME CHECKS ##
    if(multicore && !require(multicore)) stop("multicore package requested but not installed")
    if(multicore && is.null(n.cores)){
        n.cores <- parallel:::detectCores()
    }


    ## COMPUTATIONS ##
    if(multicore){
        res <- mclapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, w.dens=w.dens, w.trunc=w.trunc,
                                                                 init.tree=init.tree, n.iter=n.iter, sample.every=sample.every,
                                                                 tune.every=tune.every, pi.param1=pi.param1, pi.param2=pi.param2,
                                                                 quiet=TRUE),
                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
    } else {
         res <- lapply(1:n.runs, function(i)  outbreaker(dna=dna, dates=dates, w.dens=w.dens, w.trunc=w.trunc,
                                                                 init.tree=init.tree, n.iter=n.iter, sample.every=sample.every,
                                                                 tune.every=tune.every, pi.param1=pi.param1, pi.param2=pi.param2,
                                                         quiet=TRUE))
    }


    ## NAME RUNS ##
    names(res) <- paste("run",1:n.runs,paste="")


    ## RETURN RESULTS ##
    return(res)
} # end outbreaker.parallel

