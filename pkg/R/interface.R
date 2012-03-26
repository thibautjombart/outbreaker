#############
## outbreaker
#############
outbreaker <- function(
                      file.sizes="out-popsize.txt", file.sample="out-sample.txt"){

    

}




testrawtoC <- function(x){
    if(!require(ape)) stop("ape package is required.")

    n <- as.integer(nrow(x))
    p <- as.integer(ncol(x))
    x <- unlist(as.list(x),use.names=FALSE)

    .C("DNAbin2list_dnaseq", x, n, p, PACKAGE="outbreaker")
    return()
}
