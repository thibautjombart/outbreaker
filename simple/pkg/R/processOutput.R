get.TTree.simple <- function(x, burnin=1e5){
    if(all(x$chains$step<=burnin)) stop("requested burn-in exeeds the number of chains")

    ## CREATE OUTPUT LIST ##
    res <- list()
    res$idx <- 1:length(x$collec.dates)
    res$collec.dates <- x$collec.dates


    ## PROCESS CHAINS ##
    chains <- x$chains[x$chains$step>burnin, ] # discard burn-in

    ## get ancestors ##
    temp <- chains[,grep("alpha",names(chains))]
    res$ances <- apply(temp,2, function(e) names(table(e))[which.max(as.numeric(table(e)))])
    res$ances <- as.integer(res$ances)
    res$ances[res$ances<1] <- NA

    ## get infection dates ##
    temp <- chains[,grep("Tinf",names(chains))]
    res$inf.dates <- apply(temp,2,median)

    ## get probabilities of ancestries ##
    temp <- chains[,grep("alpha",names(chains))]
    res$p.ances <- apply(temp,2, function(e) max(as.numeric(table(e)))/sum(table(e)))

    ## get ancestor->descendent mutations ##
    D <- as.matrix(x$D)
    res$nb.mut <- sapply(1:length(res$idx), function(i) D[res$idx[i],res$ances[i]])

    ## get kappa ##
    temp <- chains[,grep("kappa",names(chains))]
    res$n.gen <- apply(temp,2, function(e) names(table(e))[which.max(as.numeric(table(e)))])
    res$n.gen <- as.integer(res$n.gen)
    res$p.gen <- apply(temp,2, function(e) max(as.numeric(table(e)))/sum(table(e)))

    ## get infectiousness curves ##
    timeSpan <- c(min(res$inf.dates),max(res$inf.dates)+length(x$w))
    f1 <- function(infDate, w){
        dens <- rep(0,diff(timeSpan)+1)
        idxStart <- infDate-timeSpan[1]+1
        dens[idxStart:(idxStart+length(w)-1)] <- w
        dates <- seq(timeSpan[1],timeSpan[2])
        res <- data.frame(dates,dens)
        return(as.matrix(res))
    }

    res$inf.curves <- lapply(res$inf.dates, f1, x$w)


    ## SET CLASS AND RETURN ##
    class(res) <- "TTree.simple"
    return(res)
} # end get.TTree.simple






#############
## as.igraph
#############
as.igraph.TTree.simple <- function(x, arr.width=c("proba", "mutations"), ...){
    if(!require(igraph)) stop("package igraph is required for this operation")
    arr.width <- match.arg(arr.width)

     ## GET DAG ##
    from <- x$ances
    to <- x$idx
    isNotNA <- !is.na(from) & !is.na(to)
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    vnames <- as.character(unique(unlist(dat)))
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=x$inf.dates[vnames]))

    ## SET WEIGHTS ##
    E(out)$weight <- x$nb.mut[isNotNA]

    ## SET ARROW WIDTH ##
    if(arr.width=="proba"){
        E(out)$width <- round(4*x$p.ances)+1
    }

    if(arr.width=="mutations"){
        temp <- max(E(out)$weight) - E(out)$weight
        temp <- temp/max(temp) * 4
        E(out)$width <- round(temp)+1
    }
    
    return(out)
}
