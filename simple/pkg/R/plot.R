#####################
## plot.TTree.simple
#####################
plot.TTree.simple <- function(x, y=NULL, edge.col="black", col.edge.by="prob",
                              col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/", ...){
    if(!require(igraph)) stop("igraph is required")
    if(!require(adegenet)) stop("adegenet is required")
    if(!inherits(x,"TTree.simple")) stop("x is not a TTree.simple object")
    if(!col.edge.by %in% c("dist","n.gen","prob")) stop("unknown col.edge.by specified")

    ## get graph ##
    g <- as.igraph(x, edge.col=edge.col, col.edge.by=col.edge.by, col.pal=col.pal, annot=annot, sep=sep)

     ## make plot ##
    plot(g, layout=attr(g,"layout"), ...)

    ## return graph invisibly ##
    return(invisible(g))

} # end plot.TTree.simple







#############
## epicurves
#############
epicurves <- function (x, col=NULL, bg="lightgrey", line.col="white", coef=1, max.arr=5,...) {
    if(!require(adegenet)) stop("adegenet is required")
    if(!inherits(x,"TTree.simple")) stop("x is not a TTree.simple object")

    ## GET USEFUL INFO ##
    N <- length(x$idx)
    timeSpan <- range(x$inf.curves[[1]][,1])
    if(is.null(col)){
        colPal <- colorRampPalette(c("grey30","blue","green3","orange","red","brown4","purple"))
        col <- colPal(N)
    }


    ## MAKE EMPTY PLOT ##
    plot(0,0,type="n",xlim=timeSpan+c(-1,1), ylim=c(0,N+1), xlab="Dates",ylab="Individual index",...)
    rect(timeSpan[1]-2,-2,timeSpan[2]+2,N+2, col=bg)
    abline(h=1:N, col="white",lwd=3)
    abline(h=1:N, col=transp(col),lwd=2)
    abline(v=pretty(timeSpan,4), col=line.col)

    ## DRAW INFECTIOUSNESS CURVES ##
    for(i in 1:N){
        temp <- x$inf.curves[[i]][x$inf.curves[[i]][,2]> 1e-12,,drop=FALSE]
        x.coord <- c(temp[,1], rev(temp[,1]))
        y.coord <- c(i+temp[,2]*coef, rep(i,nrow(temp)))
        polygon(x.coord, y.coord, col=transp(col[i]),border=col[i])
        points(temp[,1], i+(temp[,2]*coef), type="o", pch=20,cex=0.5, col=col[i])
    }


    ## ADD COLLECTION DATES ##
    points(x$collec.dates, 1:N, pch=20, cex=2, col=col)


    ## ADD INFECTIONS DATES ##
    points(x$inf.dates, 1:N, cex=2, col=col)


    ## ADD INFECTIONS ##
    arr.w <- (x$p.ances- 1/(N-1)) * max.arr
    arr.w[arr.w<0.5] <- 0.5
    arrows(x$inf.date,x$ances, x$inf.dates, x$idx, angle=15, col=col[x$ances], lwd=arr.w)

} # end epicurves











##############
## plot.chains
##############
plot.chains <- function(x, what="post", type=c("series","density"), omit.first=0, dens.all=TRUE,
                        col=rainbow(x$n.runs), lty=1, lwd=1, main=what, ...){
    ## HANDLE ARGUMENTS ##
    type <- match.arg(type)
    n.runs <- x$n.runs
    if(!what %in% names(x$chains)) stop(paste(what,"is not a column of x$chains"))
    if(!is.null(col)) col <- rep(col, length = n.runs)
    if(!is.null(lty)) lty <- rep(lty, length = n.runs)
    if(!is.null(lwd)) lwd <- rep(lwd, length = n.runs)
    if(is.null(omit.first)){
        omit.first <- max(res$burnin, res$find.import.at, res$tune.end)
    }

    ## GET DATA TO PLOT ##
    dat <- cbind(x$chains$step[x$chains$run==1],data.frame(split(x$chains[,what], x$chains$run)))
    names(dat) <- c("step", paste(what, 1:n.runs,sep=""))
    if(!any(dat$step>omit.first)) stop("omit.first is greater than the number of steps in x")
    dat <- dat[dat$step>omit.first,,drop=FALSE]

    ## MAKE PLOT ##
    if(type=="series"){
        matplot(dat$step, dat[,-1,drop=FALSE], type="l", col=col, lty=lty, xlab="MCMC iteration", ylab="value", main=main, ...)
    }

    if(type=="density"){
        ## add general density if needed ##
        temp <- lapply(dat[, -1, drop=FALSE], density)
        if(dens.all){
            temp[[n.runs+1]] <- density(unlist(dat[,-1,drop=FALSE]))
            col <- c(col, "black")
            lty <- c(lty, 1)
            lwd <- c(lwd, 3)
            n.runs <- n.runs+1
        }
        range.x <- range(sapply(temp, function(e) e$x))
        range.y <- range(sapply(temp, function(e) e$y))
        plot(1,type="n", xlim=range.x, ylim=range.y, xlab="value", ylab="density", main=main, ...)
        invisible(lapply(1:n.runs, function(i) lines(temp[[i]], col=col[i], lty=lty[i], lwd=lwd[i])))
    }

    return(invisible())
} # end plot.chains






##############
## transGraph
##############
transGraph <- function(x, labels=NULL, omit.first=x$burnin, threshold=0.2, col.pal=NULL, curved.edges=TRUE, ...){
    ## CHECKS ##
    if(!require(igraph)) stop("igraph is required")
    if(!require(adegenet)) stop("adegenet is required")

    ## HANDLE ARGUMENTS ##
    if(omit.first> max(x$chains$step)) stop("omit.first exceeds the number of chains in the output")
    if(is.null(col.pal)){
        col.pal <- function(n){
            return(grey(seq(1,0,length=n)))
        }
    }

    ## GET ANCESTRY DATA ##
    ances <- x$chains[x$chains$step>=omit.first, grep("alpha",names(x$chains)),drop=FALSE]
    tabances <- apply(ances,2,table)
    N <- ncol(ances)

    ## GET DATA FOR GRAPH ##
    ## ancestor, descendents, list of nodes
    to.old <- rep(1:N, sapply(tabances,length))
    from.old <- as.numeric(unlist(lapply(tabances,names)))
    from.old[from.old<1] <- NA
    to.old[to.old<1] <- NA
    isNotNA <- !is.na(from.old) & !is.na(to.old)
    vnames <- sort(unique(c(from.old,to.old)))
    from <- match(from.old,vnames)
    to <- match(to.old,vnames)

    ## support for the ancestries
    support <- unlist(lapply(tabances, function(e) e/sum(e)))[isNotNA]
    edge.col <- num2col(support, col.pal=col.pal, x.min=0, x.max=1)

    ## average dates of infection
    Tinf <- x$chains[x$chains$step>=omit.first, grep("Tinf",names(x$chains)),drop=FALSE]
    inf.dates <- apply(Tinf,2,mean)
    names(inf.dates) <- 1:N

    ## remove weakly supported ancestries
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    if(sum(support>=threshold)==0) warning("threshold to high - no edge left")
    dat <- dat[support>=threshold, , drop=FALSE]

    ## convert to graph
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=inf.dates[as.character(vnames)]))


    ## SET PARAMETERS OF THE GRAPH ##
    ## edges
    E(out)$color <- edge.col[support>=threshold]
    E(out)$support <- support[support>=threshold]
    E(out)$curved <- curved.edges

    ## set layout
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(minx=V(out)$dates, maxx=V(out)$dates), rescale=FALSE)


    ## MAKE THE PLOT ##
    plot(out, layout=attr(out, "layout"), ...)


    ## RETURN OBJECT ##
    return(invisible(out))

} # end transGraph
