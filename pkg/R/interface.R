
#######################
## auxiliary functions
#######################

outbreaker <- function(x, duration=100, weightNaGen=0.001){
    if(!inherits(x, "outbreaker-data")) stop("x is not an outbreaker-data object.")

    ## int *nbPatients, int *duration, int *nbAdmVec, int *nbPosSwab, int *nbNegSwab, int *nbColPatients, int *idxColPatients, int *nbSeqPat, int *wardVec, int *tAdmVec, int *tDisVec, int *tPosSwab, int *tNegSwab, int *hospPres, int *idxSeqVec, double *tCollecVec, unsigned char *DNAbinInput, int *totNbSeq, int *seqLength, double *weightNaGen
    nbPatients <- as.integer(length(x$n.swab))
    duration <- as.integer(duration)
    nbAdmVec <- as.integer(unlist(x$n.adm))
    nbPosSwab <- as.integer(sapply(x$swab, function(e) sum(e==1)))
    nbNegSwab <- as.integer(sapply(x$swab, function(e) sum(e==0)))
    nbColPatients <- as.integer(sum(nbPosSwab>0))
    idxColPatients <- as.integer(which(nbPosSwab>0)-1) # -1 important!
    nbSeqPat <- as.integer(x$n.seq)
    wardVec <- as.integer(x$ward)
    tAdmVec <- as.double(unlist(x$t.adm))
    tDisVec <- as.double(unlist(x$t.dis))
    tPosSwab <- as.double(unlist(x$t.swab)[unlist(x$swab)==1])
    tNegSwab <- as.double(unlist(x$t.swab)[unlist(x$swab)==0])
    hospPres <- as.double(unlist(x$hosp.pres))
    idxSeqVec <- as.integer(unlist(x$idx.dna)-1) # -1 important!
    tCollecVec <- as.double(unlist(x$t.dna))
    DNAbinInput <- unlist(as.list(x$dna),use.names=FALSE)
    totNbSeq <- as.integer(length(x$dna))
    seqLength <- as.integer(length(x$dna[[1]]))
    weightNaGen <- as.double(weightNaGen)

    .C("R_outbreaker", nbPatients, duration, nbAdmVec, nbPosSwab, nbNegSwab, nbColPatients, idxColPatients, nbSeqPat, wardVec, tAdmVec, tDisVec, tPosSwab, tNegSwab, hospPres, idxSeqVec, tCollecVec, DNAbinInput, totNbSeq, seqLength, weightNaGen, PACKAGE="outbreaker")

    cat("\nComputation finished.")
    return(invisible())
}

#############
## outbreaker
#############
## outbreaker <- function(
##                       file.sizes="out-popsize.txt", file.sample="out-sample.txt"){



## }



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



## testSimLike <- function(mu=0.01, p.new.patient=0.5, seq.length=100, min.n=10,...){
##     x <- sim.mrsa(mu=mu, p.new.patient=p.new.patient, seq.length=seq.length, ...)
##     while(nrow(x$dna)<min.n){
##         x <- sim.mrsa(mu=mu, p.new.patient=p.new.patient, seq.length=seq.length, ...)
##     }

##     ## get combinations of patient->patient relationships
##     npatients <- max(x$ances.pat,na.rm=TRUE)
##     temp <- lapply(1:npatients, function(i) rbind(from=setdiff(1:npatients,i),to=rep(i,npatients-1)))
##     allCombi <- Reduce(cbind,temp) # all possible infections

##     res <- double(0)

##     for(i in 1:ncol(allCombi)){
##         from <- allCombi[1,i]
##         to <- allCombi[2,i]
##         si <- which(x$patient==to)
##         sj <- which(x$patient==from)
##         ti <- as.double(x$date[si])
##         tj <- as.double(x$date[sj])
##         si <- as.integer(si-1)
##         sj <- as.integer(sj-1)



##         ##test_genlike(unsigned char *DNAbinInput, int *n, int *length, int *s_i, , int *s_j, double *t_i, double *t_j, int *m_i, int *m_j, double *nu1, double *nu2, double *alpha, double *tau, double *out)

##         llik <- unlist(.C("test_genlike", unlist(as.list(x$dna),use.names=FALSE), nrow(x$dna), ncol(x$dna),
##                    si, sj, ti, tj, length(si), length(sj), as.double(mu), as.double(mu),
##                    as.double(p.new.patient), as.double(1), double(1), PACKAGE="outbreaker")[14])

##         res <- c(res,llik)
##     }

##     ## identify actual infections in allCombi ##
##     candi <- apply(allCombi,2,paste, collapse="-")
##     real <- apply(x$ances.pat,1,paste, collapse="-")


##     res <- data.frame(from=allCombi[1,],to=allCombi[2,],actual=candi %in% real,likelihood=res)
##     return(res)
## }

## ##
## ##  library(outbreaker)
## ##  test <- testSimLike()
## ##  library(ggplot2)
## ##  qplot(likelihood,data=test,colour=actual, geom="density")
## ##
## ##
