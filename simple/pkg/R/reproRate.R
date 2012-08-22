##########
## get.Rt
##########
get.Rt <- function(x, burnin=1e5, plot=TRUE, type=c("CI", "boxplot", "lines"), lines=TRUE,
                   CI.level=0.05, fill.col="gold", lines.col=transp("grey"), ...){
    if(!require(adegenet)) stop("the adegenet package is required.")
    type <- match.arg(type)
    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]

    ## get ancestries
    ances <- dat[,grep("alpha", names(x$chains)),drop=FALSE]
    tabAnces <- apply(ances, 1, table)

    ## get infection times
    Tinf <-  dat[,grep("Tinf", names(x$chains)),drop=FALSE]
    timeSpan <- range(Tinf)
    timeStep <- seq(timeSpan[1],timeSpan[2],by=1)
    emptyOut <- rep(0, length(timeStep))
    names(emptyOut) <- timeStep

    ## function to get Rt for one chain
    f1 <- function(i){
        e <- tabAnces[[i]][-1] # -1: remove '0's
        out <- emptyOut
        Tinf.temp <- Tinf[i,as.integer(names(e))]
        nbCasePerTimeStep <- tapply(as.numeric(e),as.numeric(Tinf.temp),sum)
        out[names(nbCasePerTimeStep)] <- nbCasePerTimeStep
        return(out)
    }

    ## GET RT FOR ALL RELEVANT CHAINS ##
    res <- lapply(1:nrow(dat), f1)
    res <- t(Reduce("rbind",res))


    ## MAKE PLOT IF NEEDED ##
    if(plot){
        if(type=="CI"){
            ## CI-based
            stat <- apply(res, 1, quantile, c(CI.level,0.5, 1-CI.level))
            xcoord <- as.numeric(colnames(stat))
            matplot(xcoord, t(stat), type="n", ...)
            polygon(c(xcoord,rev(xcoord)), c(stat[3,], rev(stat[1,])), col=fill.col, border=NA)
            if(lines){
                matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
            }
            matplot(xcoord, t(stat), type="l", lty=1, col="black", lwd=2, add=TRUE)
        }



        ## boxplot-based
        if(type=="boxplot"){
            boxplot(t(res), col=fill.col, at=as.integer(rownames(res)),
                    xlab="Time", ylab="Effective reproduction number (R(t))", ...)
            if(lines){
                matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
            }
        }

        ## just lines
        if(type=="lines"){
            matplot(rownames(res),res, type="l", lty=1, col=lines.col,
                    xlab="Time", ylab="Effective reproduction number (R(t))", ...)
        }
    }

    return(res)
} # end get.Rt
