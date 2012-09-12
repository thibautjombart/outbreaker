##########
## get.Rt
##########
get.Rt <- function(x, burnin=2e4, plot=TRUE, type=c("boxplot", "lines"), lines=FALSE,
                   CI.level=0.05, fill.col="gold", lines.col=transp("grey"), ...){
    if(!require(adegenet)) stop("the adegenet package is required.")
    type <- match.arg(type)

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]

    ## get ancestries
    ances <- dat[,grep("alpha", names(x$chains)),drop=FALSE] # table of ancestries
    tabAnces <- apply(ances, 1, table) # count nb of descendents per case for each chain

    ## get infection times
    Tinf <-  dat[,grep("Tinf", names(x$chains)),drop=FALSE]
    timeSpan <- range(Tinf)
    timeStep <- seq(timeSpan[1],timeSpan[2],by=1)
    emptyOut <- rep(NA, length(timeStep))
    names(emptyOut) <- timeStep

    ## function to get Rt for one chain
    f1 <- function(i){
        ## get nb of descendents per ancestor
        e <- tabAnces[[i]][-1] # -1: remove '0's

        ## create empty output
        out <- emptyOut

        ## get infection times of ancestors
        Tinf.temp <- Tinf[i,as.integer(names(e))]

        ## count mean nb of secondary infections created by cases infected at each time step
        meanNbCasePerTimeStep <- tapply(as.numeric(e),as.numeric(Tinf.temp),mean)
        out[names(meanNbCasePerTimeStep)] <- meanNbCasePerTimeStep
        return(out)
    }

    ## GET RT FOR ALL RELEVANT CHAINS ##
    res <- lapply(1:nrow(dat), function(i) f1(i))
    res <- t(Reduce("rbind",res))


    ## MAKE PLOT IF NEEDED ##
    if(plot){
        ## if(type=="CI"){
        ##     ## CI-based
        ##     stat <- apply(res, 1, quantile, c(CI.level,0.5, 1-CI.level), na.rm=TRUE)
        ##     xcoord <- as.numeric(colnames(stat))
        ##     matplot(xcoord, t(stat), type="n", ...)
        ##     polygon(c(xcoord,rev(xcoord)), c(stat[3,], rev(stat[1,])), col=fill.col, border=NA)
        ##     if(lines){
        ##         matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
        ##     }
        ##     matplot(xcoord, t(stat), type="l", lty=1, col="black", lwd=2, add=TRUE)
        ## }



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







#############
## get.incid
#############
get.incid <- function(x, burnin=2e4, plot=TRUE, type=c("CI", "boxplot", "lines"), lines=FALSE,
                   CI.level=0.05, fill.col="gold", lines.col=transp("grey"), ...){
    if(!require(adegenet)) stop("the adegenet package is required.")
    type <- match.arg(type)

    ## GET DATA ##
    ## remove burnin
    if(!any(x$chains$step>burnin)) stop("burnin too high - no chain left")
    dat <- x$chains[x$chains$step>burnin,,drop=FALSE]

    ## get infection times
    Tinf <-  as.matrix(dat[,grep("Tinf", names(x$chains)),drop=FALSE])
    timeSpan <- range(Tinf)
    timeStep <- seq(timeSpan[1],timeSpan[2],by=1)
    emptyOut <- rep(0, length(timeStep))
    names(emptyOut) <- timeStep

    ## function to get Rt for one chain
    f1 <- function(i){
        ## get nb of descendents per ancestor
        e <- Tinf[i,,drop=TRUE] # -1: remove '0's

        ## create empty output
        out <- emptyOut

        ## fill in output
        out[names(table(e))] <- table(e)

        return(out)
    }

    ## GET RT FOR ALL RELEVANT CHAINS ##
    res <- lapply(1:nrow(dat), function(i) f1(i))
    res <- t(Reduce("rbind",res))


    ## MAKE PLOT IF NEEDED ##
    if(plot){
        if(type=="CI"){
            ## CI-based
            stat <- apply(res, 1, quantile, c(CI.level,0.5, 1-CI.level), na.rm=TRUE)
            xcoord <- as.numeric(colnames(stat))
            matplot(xcoord, t(stat), type="n", ylab="Incidence", ...)
            polygon(c(xcoord,rev(xcoord)), c(stat[3,], rev(stat[1,])), col=fill.col, border=NA)
            if(lines){
                matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
            }
            matplot(xcoord, t(stat), type="l", lty=1, col="black", lwd=2, add=TRUE)
        }



        ## boxplot-based
        if(type=="boxplot"){
            boxplot(t(res), col=fill.col, at=as.integer(rownames(res)),
                    xlab="Time", ylab="Incidence", ...)
            if(lines){
                matplot(rownames(res),res, type="l", lty=1, col=lines.col, lwd=1, add=TRUE)
            }
        }

        ## just lines
        if(type=="lines"){
            matplot(rownames(res),res, type="l", lty=1, col=lines.col,
                    xlab="Time", ylab="Incidence", ...)
        }
    }

    return(res)
} # end get.incid
