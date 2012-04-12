
#######################
## auxiliary functions
#######################

## typedef struct{
##     int NbPatients, T;
##     int * NbAdmissions;
##     int * NbPosSwabs;
##     int * NbNegSwabs;
##     int NbColonisedPatients; /* those who have at least one positive swab */
##     int * indexColonisedPatients; /* those who have at least one positive swab */
##     int * M; /* number of sequences in each patient */
## } nb_data;

.prep.nbdata <- function(){
    
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
