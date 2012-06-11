
#######################
## auxiliary functions
#######################

outbreaker <- function(dna, dates){
    ## CHECKS ##
    if(!require(ape)) stop("the ape package is required but not installed")
    if(!inherits(dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if(!is.matrix(dna)) dna <- as.matrix(dna)
    if(is.character(dates)) stop("dates are characters; they must be integers or dates with POSIXct format (see ?as.POSIXct)")


    ## PROCESS INPUTS ##
    ## dna ##
    n.ind <- as.integer(nrow(dna))
    n.nucl <- as.integer(ncol(dna))
    dna <- unlist(as.list(dna),use.names=FALSE)
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

    ## make sure minimum date is 0 ##
    dates <- as.integer(dates-min(dates))

    ## .C("Rinput2data", dna, dates, n.ind, n.nucl, PACKAGE="outbreaker")
    .C("R_outbreaker", dna, dates, n.ind, n.nucl, PACKAGE="outbreaker")

    cat("\nComputations finished.")
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


