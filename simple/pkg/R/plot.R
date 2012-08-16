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
