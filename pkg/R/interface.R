#############
## outbreaker
#############
## outbreaker <- function(
##                       file.sizes="out-popsize.txt", file.sample="out-sample.txt"){



## }



##################
## test functions
##################
testrawtoC <- function(x){
    if(!require(ape)) stop("ape package is required.")

    n <- as.integer(nrow(x))
    p <- as.integer(ncol(x))
    x <- unlist(as.list(x),use.names=FALSE)

    .C("DNAbin2list_dnaseq", x, n, p, PACKAGE="outbreaker")
    return()
}



testSimLike <- function(mu=0.01, p.new.patient=0.5, seq.length=100, ...){
    x <- sim.mrsa(mu=mu, p.new.patient=p.new.patient, seq.length=seq.length, ...)
    while(nrow(x$dna)<10){
        x <- sim.mrsa(mu=mu, p.new.patient=p.new.patient, seq.length=seq.length, ...)
    }

    from <- 1L
    to <- 2L
    si <- which(x$patient==to)
    sj <- which(x$patient==from)
    ti <- as.double(x$date[si])
    tj <- as.double(x$date[sj])
    res <- double(1)
    si <- as.integer(si-1)
    sj <- as.integer(sj-1)


##test_genlike(unsigned char *DNAbinInput, int *n, int *length, int *s_i, , int *s_j, double *t_i, double *t_j, int *m_i, int *m_j, double *nu1, double *nu2, double *alpha, double *tau, double *out)

    temp <- .C("test_genlike", unlist(as.list(x$dna),use.names=FALSE), nrow(x$dna), ncol(x$dna),
               si, sj, ti, tj, length(si), length(sj), as.double(mu), as.double(mu),
               as.double(p.new.patient), as.double(1), double(1), PACKAGE="outbreaker")
    res <- temp[14]
}

##
##  library(outbreaker)
##
##
##
