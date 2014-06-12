


###########
## simMers
###########
##
## R0: basic repro number
## infec.curve: generation time distribution
## disp: mean distance for host movement
simMers <- function(R0, n.patches=1, xy.patches=NULL,
                    d.gentime, d.incubation=d.gentime, disp.patches=NULL,
                    rate.transfer.case=0, rate.import.case=0, mean.dispersal=0, names=NULL,
                    max.time=100, max.infected=NULL, grid.size=10){

    ## HANDLE ARGUMENTS ##
    if(!is.null(xy.patches)){
        if(is.data.frame(xy.patches)) xy.patches <- as.matrix(xy.patches)
        if(!is.matrix(xy.patches)) stop("xy.patches is not a matrix")
        if(nrow(xy.patches)!=n.patches) stop("xy.patches must have one row per patch")
    }
    if(is.null(names)) names <- paste("patch",1:n.patches, sep="_")

    ## set max.time
    if(is.null(max.time)) {
        max.time <- 5e3
    }

    ## set max.infected
    if(is.null(max.infected)) {
        max.infected <- 1e4
    }

    ## spatial coordinates of the patches
    if(is.null(xy.patches)) xy.patches <- matrix(runif(2*n.patches, min=0, max=grid.size), ncol=2)
    colnames(xy.patches) <- c("x","y")

    ## distances between patches
    if(is.null(disp.patches)) disp.patches <- as.matrix(dist(xy.patches))
    disp.patches <- as.matrix(disp.patches)

    ## scale distributions
    d.gentime <- d.gentime/sum(d.gentime)
    d.incubation <- d.incubation/sum(d.incubation)

    ## add long tails to distributions
    d.gentime <- c(d.gentime, rep(0, max.time))
    d.incubation <- c(d.incubation, rep(0, max.time))



    ## AUXILIARY FUNCTIONS ##
    ## function to build one patch ##
    makePatch <- function(xy=NULL){
        ## find spatial coords if not provided
        out <- list(Ninf=1, Tinf=0, xy=xy, type="local")
        return(out)
    }

    ## add cases in a patch ##
    addCases <- function(patch, n, type="local"){
        ## add dates of infection for new infected
        patch$Tinf <- c(patch$Tinf, rep(t, n))

        ## add new infected to counter
        patch$Ninf <- patch$Ninf + n

        ## add new type
        patch$type <- c(patch$type, rep(type, n))

        ## return
        return(patch)
    }

    ## compute global force of infection in a patch
    ## t is current time
    forceOfInfection <- function(patch, t){
        out <- sum(R0*d.gentime[1+t-patch$Tinf])
        return(out)
    }

    ## function defining new infections in a patch ##
    ## t is the time
    spreadInPatch <- function(patch, t){
        ## global force of infection (FOI)
        global.FOI <-forceOfInfection(patch, t)

        ## draw number of new infected
        new.inf <- rpois(1, lambda=global.FOI)

        ## add dates of infection for new infected
        patch <- addCases(patch, new.inf, type="local")

        ## return
        return(patch)
    }

    ## migration between patches ##
    migration <- function(patches, t){
        ## escape if just one patch
        if(length(patches)==1) return(patches)

        ## compute force of dispersal
        ## (Ninf * disp.rate, for each patch)
        patches.FOI <- sapply(patches, forceOfInfection, t)
        disp.force <- patches.FOI * disp.patches * rate.transfer.case

        ## matrix of nb of trans-hospital infections
        ## row = from, column = to
        nb.trans.inf <- matrix(sapply(as.vector(disp.force), rpois, n=1),ncol=n.patches)

        ## nb of new cases in each patch
        nb.newcases <- apply(nb.trans.inf,2,sum)

        ## add cases to the patches
        patches <- lapply(1:n.patches, function(i) addCases(patches[[i]], nb.newcases[i], type="transfer"))

        ## return
        return(patches)
    }

    ## importation from reservoir ##
    importation <- function(patches){
        ## number of new imported cases per patch
        nb.import <- rpois(n.patches, rate.import.case)

        ## add cases to the patches
        patches <- lapply(1:n.patches, function(i) addCases(patches[[i]], nb.newcases[i], type="imported"))

        ## return
        return(patches)
    }



    ## INITIALIZE PATCHES ##
    patches <- lapply(1:n.patches, function(i) makePatch(xy.patches[i,]))
    names(patches) <- names


    ## RUN SIMULATIONS ##
    ## continue counter
    CONTINUE <- TRUE

    ## time
    t <- 0

    ## output: incidence of local cases
    incid.local <- n.patches

    ## output: incidence of cross-patch infections
    incid.transfer <- 0

    ## output: incidence of cross-patch infections
    incid.import <- 0

    ## run simulations
    while(CONTINUE){
        ## this is a beautiful day, full of opportunities...
        t <- t+1

        ## spread disease within patches
        patches <- lapply(patches, spreadInPatch, t)

        ## trans-patches infections
        patches <- migration(patches, t)

        ## compute total number of infected
        total.infected <- sum(sapply(patches, function(e) e$Ninf))

        ## compute local incidence
        incid.local <- c(incid.local, sum(sapply(patches, function(e) sum(e$Tinf==t & e$type=="local"))))

        ## compute 'transfer' (across patch) incidence
        incid.transfer <- c(incid.transfer, sum(sapply(patches, function(e) sum(e$Tinf==t & e$type=="transfer"))))

        ## compute 'imports' (from reservoir) incidence
        incid.import <- c(incid.import, sum(sapply(patches, function(e) sum(e$Tinf==t & e$type=="imported"))))

        ## compute
        ## update 'CONTINUE'
        if(t>=max.time ||total.infected>=max.infected) CONTINUE <- FALSE
    }


    ## SHAPE AND RETURN OUTPUT ##
    incidence <- data.frame(local=incid.local,
                            transfer=incid.transfer,
                            imported=incid.import)

    res <- list(incidence=incidence,
                nb.cases=total.infected,
                xy=xy.patches,
                patches=patches)


    return(res)
} # end simMers
