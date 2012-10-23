################################################
## Function which looks for simulation results
## recursively in all directories and compiles
## results in a file

poolResults <- function(dir=getwd(), file="pooledResults.csv"){

    ## get input files ##
    in.files <- list.files(path=dir, include=TRUE, pattern="in.csv", recursive=TRUE)


    ## get output files ##
    out.files <- list.files(path=dir, include=TRUE, pattern="out.csv", recursive=TRUE)


    ## keep only files with inputs and outputs ##
    temp1 <- sub("in.csv", "", in.files)
    temp2 <- sub("out.csv", "", out.files)
    temp <- table(c(temp1,temp2))
    toKeep <- names(temp)[temp==2]
    in.files <- in.files[in.files %in% paste(toKeep,"in.csv",sep="")]
    out.files <- out.files[out.files %in% paste(toKeep,"out.csv",sep="")]

    ## compile results ##
    res <- list()
    for(i in 1:length(in.files)){
        temp <- cbind.data.frame(read.csv(in.files[i], stringsAsFactors=FALSE),
                                 read.csv(out.files[i], stringsAsFactors=FALSE),
                             stringsAsFactors=FALSE)
        temp$X <- NULL
        res <- rbind.data.frame(res, temp)
    }

    ## write res to file ##
    write.csv(res, file=file)

    ## return res ##
    return(res)
} # end poolResults
