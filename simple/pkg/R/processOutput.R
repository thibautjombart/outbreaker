####################
## get.TTree.simple
####################
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
as.igraph.TTree.simple <- function(x, edge.col="black", col.edge.by="prob",
                              col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/", ...){
    if(!require(igraph)) stop("package igraph is required for this operation")
    if(!require(adegenet)) stop("adegenet is required")
    if(!inherits(x,"TTree.simple")) stop("x is not a TTree.simple object")
    if(!col.edge.by %in% c("dist","n.gen","prob")) stop("unknown col.edge.by specified")

    ## GET DAG ##
    from.old <- x$ances
    to.old <- x$id
    isNotNA <- !is.na(from.old) & !is.na(to.old)
    vnames <- sort(unique(c(from.old,to.old)))
    from <- match(from.old,vnames)
    to <- match(to.old,vnames)
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=x$inf.dates[vnames]))

    ## SET VARIOUS INFO ##
    E(out)$dist <- x$nb.mut[isNotNA]
    E(out)$prob <- x$p.ances[isNotNA]
    E(out)$n.gen <- x$n.gen[isNotNA]
    E(out)$p.kappa <- x$p.gen[isNotNA]

   ## SET EDGE COLORS ##
    if(is.null(col.pal)){
        col.pal <- function(n){
            return(grey(seq(0.75,0,length=n)))
        }
    }
    if(col.edge.by=="prob") edge.col <- num2col(E(out)$prob, col.pal=col.pal, x.min=0, x.max=1)
    if(col.edge.by=="dist") edge.col <- num2col(E(out)$dist, col.pal=col.pal, x.min=0, x.max=1)
    if(col.edge.by=="n.gen") edge.col <- num2col(E(out)$n.gen, col.pal=col.pal, x.min=0, x.max=1)

    E(out)$color <- edge.col

    ## SET EDGE LABELS ##
    n.annot <- sum(annot %in% c("dist","n.gen","prob"))
    lab <- ""
    if(!is.null(annot) && n.annot>0){
        if("dist" %in% annot) lab <- E(out)$dist
        if("n.gen" %in% annot) lab <- paste(lab, E(out)$n.gen, sep=sep)
        if("prob" %in% annot) lab <- paste(lab, round(E(out)$prob,2), sep=sep)
    }
    lab <- sub(paste("^",sep,sep=""),"",lab)
    E(out)$label <- lab


    ## SET LAYOUT ##
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(minx=x$inf.dates, maxx=x$inf.dates))

    return(out)
} # end as.igraph.TTree.simple
