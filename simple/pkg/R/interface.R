
#######################
## auxiliary functions
#######################

outbreaker <- function(dna, dates, w.type=1, w.param=c(2,0,0), w.trunc=15,
                       init.tree=c("seqTrack","random"),
                       n.iter=1e5, sample.every=1000, tune.every=1000, quiet=FALSE){
    ## CHECKS ##
    if(!require(ape)) stop("the ape package is required but not installed")
    if(!inherits(dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if(!is.matrix(dna)) dna <- as.matrix(dna)
    if(is.character(dates)) stop("dates are characters; they must be integers or dates with POSIXct format (see ?as.POSIXct)")
    init.tree <- match.arg(init.tree)


    ## PROCESS INPUTS ##
    ## dna ##
    n.ind <- as.integer(nrow(dna))
    n.nucl <- as.integer(ncol(dna))
    dnaraw <- unlist(as.list(dna),use.names=FALSE)
    if(n.ind != length(dates)) stop(paste("dna and dates have different number of individuals -",n.ind,"versus",length(dates)))

    ## dates in numeric ##
    if(is.numeric(dates)){
        dates <- as.integer(round(dates))
        warning("dates have been rounded to nearest integers")
    }

    ## if dates are in POSIXct format ##
    if(inherits(dates, "POSIXct")){
        dates <- as.integer(difftime(dates, min(dates), units="days"))
    }

    ## make sure minimugit pum date is 0 ##
    dates <- as.integer(dates-min(dates))

    ## parameters of generation time function ##
    w.type <- as.integer(w.type)
    w.param <- rep(w.param, length=3)
    w.param <- as.double(w.param)
    w.trunc <- as.integer(w.trunc)


    ## find initial tree ##

    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    canBeAnces <- outer(dates,dates,FUN="<=")
    diag(canBeAnces) <- FALSE

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

    ## coerce type for remaining arguments ##
    n.iter <- as.integer(n.iter)
    sample.every <- as.integer(sample.every)
    tune.every <- as.integer(tune.every)
    quiet <- as.integer(quiet)

    ## .C("Rinput2data", dna, dates, n.ind, n.nucl, PACKAGE="outbreaker") int *wType, int *wParam1, int *wParam2, int *wParam3, int *wTrunc
    .C("R_outbreaker",
       dnaraw, dates, n.ind, n.nucl,
       w.type, w.param[1], w.param[2], w.param[3], w.trunc,
       ances, n.iter, sample.every, tune.every, quiet,
       PACKAGE="outbreaker")

    cat("\nComputations finished.\n")
    return(invisible())
}



## ##################
## ## test functions
## ##################
## testrawtoC <- function(x){
##     if(!require(ape)) stop("ape package is required.")

##     n <- as.integer(nrow(x))
##     p <- as.integer(ncol(x))
##     x <- unlist(as.list(x),use.names=FALSE)

##     .C("DNAbin2list_dnaseq", x, n, p, PACKAGE="outbreaker")
##     return()
## }


