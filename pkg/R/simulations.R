####################
## sim.testrun.data
####################
##
## simple function to simulate data to test R/C interface
sim.testrun.data <- function(n.patients=10, min.n.swab=0, max.n.swab=10,
                             min.n.seq=0, max.n.seq=10, time.span=100,
                             seq.length=30){

    ## CHECKS ##
    if(!require(ape)) stop("The ape package is required.")
    if(max.n.seq>max.n.swab) max.n.seq <- max.n.swab


    ## EPI DATA ##
    ## nb of swab per patient
    n.swab <- as.integer(round(runif(n.patients, min.n.swab, max.n.swab)))

     ## simulate swab data
    swab <- lapply(n.swab, function(i) as.integer(sample(0:1,i,replace=TRUE)))
    t.swab <- lapply(n.swab, function(i) as.double(round(runif(i, 0, time.span))))

    ## simulate wards
    ward <- sapply(1:n.patients, function(i) as.integer(sample(0:1,1)))

    ## simulate admissions / discharges
    n.adm <- rep(1L, n.patients)

    t.adm <- sapply(t.swab, function(e) ifelse(all(is.na(e)), 0, min(e)-1) )
    t.adm[t.adm<0] <- 0
    t.adm <- as.double(t.adm)

    t.dis <- sapply(t.swab, function(e) ifelse(all(is.na(e)), time.span, max(e)+1) )
    t.dis[t.dis>time.span] <- time.span
    t.dis <- as.double(t.dis)

    ## presence in hospital (0:no, 1:yes)
    hosp.pres <- lapply(1:n.patients, function(i) sapply(0:time.span, function(t) ifelse(t<t.adm[i] | t>t.dis[i],0,1)))


    ## GENETIC DATA ##
    ## generate sequence from scratch
    NUCL <- c("a","t","c","g","-")
    NA.prop <- .02
    PROB <- c(rep((1-NA.prop)/4,4),NA.prop)
    seq.gen <- function(){
        res <- sample(NUCL, size=seq.length, replace=TRUE, prob=PROB)
        return(res)
    }

    ## nb of sequences per patient
    n.seq <- sapply(n.swab, function(e) round(runif(1, 0, e)))
    n.seq[n.seq<min.n.seq] <- min.n.seq
    n.seq[n.seq>max.n.seq] <- max.n.seq

    ## simulate sequences
    dna <- as.DNAbin(lapply(1:sum(n.seq), function(i) seq.gen()))
    t.dna <- mapply(function(a,b) sample(a, b, replace=TRUE), t.swab, n.seq)

    ## simulate sequence indices for the patients
    temp <- unlist(lapply(1:length(n.seq), function(i) rep(i, n.seq[i])))
    temp <- sample(temp)
    idx.dna <- lapply(1:n.patients, function(i) as.integer(which(temp==i)))


    ## FORM RESULT AND RETURN ##
    res <- list(n.swab=n.swab, swab=swab, t.swab=t.swab, ward=ward, n.adm=n.adm, t.adm=t.adm, t.dis=t.dis,
                hosp.pres=hosp.pres, n.seq=n.seq, dna=dna, t.dna=t.dna, idx.dna=idx.dna)
    class(res) <- "outbreaker-data"
    return(res)

} # end sim.testrun.data




############
## sim.mrsa
############
sim.mrsa <- function(t.max = 10, seq.length=1e4, mu=1e-5, p.new.patient=0.8,
                          nb.lineages=1, lineage.diff=50,
                          gen.time=function(){round(rexp(1,0.2))},
                          R0=function(){round(rpois(1,1.5))}){
    ## CHECKS ##
    if(!require(ape)) stop("The ape package is required.")


    ## HANDLE ARGUMENTS ##
    if(is.numeric(gen.time)){
        gen.time.val <- gen.time
        gen.time <- function(){gen.time.val}
    }

    if(is.numeric(R0)){
        R0.val <- R0
        R0 <- function(){R0.val}
    }


    ## AUXILIARY FUNCTIONS ##
    NUCL <- c("a","t","c","g")

    ## generate sequence from scratch
    seq.gen <- function(){
        res <- sample(NUCL, size=seq.length, replace=TRUE)
        return(res)
    }

    ## create substitutions for defined SNPs
    substi <- function(snp){
        res <- sapply(snp, function(e) sample(setdiff(NUCL,e),1)) # ! sapply does not work on DNAbin vectors directly
        return(res)
    }

    ## duplicate a sequence (including possible mutations)
    seq.dupli <- function(seq, T){ # T is the number of time units between ancestor and decendent
        temp <- rpois(1, seq.length*T*mu) # total number of mutations
        toChange <- sample(1:seq.length, size=temp, replace=FALSE)
        if(sum(toChange)>0) {
            seq[toChange] <- substi(seq[toChange])
        }
        return(seq)
    }

    ## duplicate a sequence with specified number of mutations
    seq.diverg <- function(seq, nmut){ # T is the number of time units between ancestor and decendent
        toChange <- sample(1:seq.length, size=nmut, replace=FALSE)
        if(sum(toChange)>0) {
            seq[toChange] <- substi(seq[toChange])
        }
        return(seq)
    }

    ## how much time before duplication occurs ?
    time.dupli <- function(){
        res <- gen.time()
        res[res<0] <- 0 # safeguard
        return(res)
    }

    ## when duplication occurs?
    date.dupli <- function(curDate){
        res <- curDate + time.dupli()
        return(res)
    }

    ## how many duplication/transmission occur?
    nb.desc <- function(){
        res <- R0()
        res[res<0] <- 0 # safeguard
        return(res)
    }

    ## what are the IDs of the (new) patients?
    patient.id <- function(cur.id, n){
        areNewPatients <- p.new.patient > runif(n,0,1)
        out <- rep(cur.id,n)
        newPatients <- seq(from=max(res$patient)+1, by=1, length=sum(areNewPatients))
        out[areNewPatients] <- newPatients
        return(out)
    }


    ## MAIN SUB-FUNCTION: REPLICATION OF ONE ISOLATE ##
    replicate.isolate <- function(seq, date, idx){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.desc()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- sapply(1:nbDes, function(i) date.dupli(date)) # find dates for descendants
        newDates <- newDates[newDates <= t.max] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq, newDates[i]-date)) # generate new sequences
        newSeq <- matrix(unlist(newSeq), byrow=TRUE, nrow=nbDes)
        res$dna <<- rbind(res$dna, newSeq) # append to general output
        res$date <<- c(res$date, newDates) # append to general output
        res$ances <<- c(res$ances, rep(idx, nbDes)) # append to general output
        res$patient <<- c(res$patient, patient.id(res$patient[idx], nbDes))
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
        return(NULL)
    }



    ## PERFORM SIMULATIONS ##
    ## INITIALIZATION
    res <- list()
    firstseq <- seq.gen()
    if(nb.lineages>1){
        temp <- lapply(nb.lineages-1, function(i) seq.diverg(firstseq,lineage.diff))
        res$dna <- matrix(c(firstseq,unlist(temp)), byrow=TRUE, nrow=nb.lineages)
    } else {
        res$dna <- matrix(firstseq, byrow=TRUE, nrow=nb.lineages)
    }
    res$patient <- 1:nb.lineages
    res$date[1:nb.lineages] <- rep(0,nb.lineages)
    res$ances[1:nb.lineages] <- rep(NA,nb.lineages)
    toExpand <- rep(TRUE,nb.lineages)

    ## SIMULATIONS: isn't simplicity beautiful?
    while(any(toExpand)){
        idx <- min(which(toExpand))
        replicate.isolate(res$dna[idx,], res$date[idx], idx)
        ##resize.result()
    }


    ## SHAPE AND RETURN OUTPUT ##
    res$ances.iso <- data.frame(from=res$ances, to=1:length(res$ances))
    res$patient <- res$patient
    res$ances.pat <- data.frame(from=res$patient[res$ances.iso$from], to=res$patient)
    res$dna <- as.DNAbin(res$dna)
    res$ances <- NULL
    D <- as.matrix(dist.dna(res$dna, model="raw")*ncol(res$dna))
    if(nrow(res$dna)>1 && mean(is.na(res$ances.iso$from))>0){
        res$mut <- mapply(function(i,j) D[i,j], res$ances.iso$from, res$ances.iso$to)
    } else {
        res$mut <- rep(NA, sum(is.na(res$ances.iso$from)))
    }
    res$call <- match.call()
    class(res) <- "simMrsa"
    return(res)
} # end sim.mrsa







############
## plotTree
############
plotTree <- function(x, ...){
    ## CHECKS ##
    if(!require(ape)) stop("The ape package is required.")
    if(!inherits(x, "simMrsa")) stop("x is not a simMrsa object")

    ## GET TREE ##
    if(nrow(x$dna)==1){
        cat("\nOnly one isolate.\n")
        return(invisible())
    }
    D <- dist.dna(x$dna, model="raw")*ncol(x$dna)
    if(all(D<1e-14)) {
        cat("\nAll isolates are identical.\n")
        return(invisible())
    }
    tre <- nj(D)

    ## MAKE THE PLOT ##
    myCol <- rainbow(max(x$patient))[x$patient]
    plot(tre, tip.col=myCol, ...)
    axisPhylo()

    return(invisible())
} # end plotTree








#############
## plotTrans
#############
plotTrans <- function(x, div.dim=3, cex=3, ...){
    ## CHECKS ##
    if(!require(ape)) stop("The ape package is required.")
    if(!require(adegenet)) stop("The adegenet package is required.")
    if(!inherits(x, "simMrsa")) stop("x is not a simMrsa object")
    if(nrow(x$dna)==1){
        cat("\nOnly one isolate.\n")
        return(invisible())
    }


    ## MAKE THE PLOT ##
    ## get colors info ##
    D <- dist.dna(x$dna, model="raw")*ncol(x$dna)
    D <- suppressWarnings(cailliez(D))
    pco1 <- dudi.pco(D,scannf=FALSE,nf=div.dim)

    ## main plot ##
    ## empty plot
    plot(x$date,x$patient, xlab="Date", ylab="Patient", type="n",xaxt="n")

    ## density
    dens <- density(x$date, bw=0.4)
    pol.x <- dens$x
    pol.y <- 0.8*max(x$patient) *dens$y/max(dens$y)
    polygon(c(pol.x,rev(pol.x)), c(pol.y,rep(0,length(pol.x))), col="grey", border="black")

    ## colorplot
    colorplot(data.frame(x$date,x$patient), pco1$li,cex=cex, defaultLevel=0.5, alpha=1, add=TRUE,...)

    ## arrows
    areNA <- is.na(x$ances.iso$from)
    x.from <- x$date[x$ances.iso$from][!areNA]
    y.from <- x$ances.pat$from[!areNA]
    x.to <- x$date[!areNA]
    y.to <- x$ances.pat$to[!areNA]
    myPal <- colorRampPalette(c("blue","red"))
    myCol <- num2col(x$mut[!areNA], col.pal=myPal)
    suppressWarnings(arrows(x.from, y.from, x.to, y.to, col=myCol, length=0.15, angle=20))

    ## legend
    val <- pretty(x$mut, n=4)
    legend("topleft", fill=num2col(val, col.pal=myPal, x.max=max(x$mut,na.rm=TRUE)), legend=val, title="nb of mutations")

    ## axis
    axis(side=1)

    return(invisible())
} # end plotTrans



