
#######################
## auxiliary functions
#######################

outbreaker <- function(dna, dates, w.dens, w.trunc=length(w.dens),
                       init.tree=c("seqTrack","random","star"),
                       n.iter=1e5, sample.every=1000, tune.every=1000,
                       pi.param1=6, pi.param2=1, quiet=FALSE){
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
    }

    ## coerce type for remaining arguments ##
    n.iter <- as.integer(n.iter)
    sample.every <- as.integer(sample.every)
    tune.every <- as.integer(tune.every)
    pi.param1 <- as.double(pi.param1)
    pi.param2 <- as.double(pi.param2)
    quiet <- as.integer(quiet)

    .C("R_outbreaker",
       dnaraw, dates, n.ind, n.nucl,
       w.dens, w.trunc,
       ances, n.iter, sample.every, tune.every, pi.param1, pi.param2, quiet,
       PACKAGE="outbreaker")

    cat("\nComputations finished.\n\n")
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


