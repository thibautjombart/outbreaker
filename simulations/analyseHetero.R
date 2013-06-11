#########################
## 2-groups simulations
#########################

## read simulation outputs
files.2g <- dir(recursive=TRUE, pattern="test2groups.csv")
dat.2g <- lapply(files.2g, read.csv)
dat.2g <- Reduce("rbind.data.frame",dat.2g)

## find true reproduction numbers
dat.files <- dir("sim2groups", recursive=TRUE, patter=".RData",full=TRUE)

findTrueR <- function(file){
    ## load file
    load(file)

    ## get ancestries
    ances <- dat$ances
    N <- length(ances)
    R <- sapply(1:N, function(i) sum(ances==i,na.rm=TRUE))
    out <- tapply(R, dat$group, mean, na.rm=TRUE)

    ## return
    return(out)
}

trueR <- t(sapply(dat.files, findTrueR))

## boxplot
boxplot(dat.2g[,c(2,3,5,6)])




#########################
## 2-groups simulations
#########################
files.zup <- dir(recursive=TRUE, pattern="zuper.txt")
dat.zup <- lapply(files.zup,read.table)
dat.zup <- Reduce("rbind.data.frame",dat.zup[sapply(dat.zup, ncol)==2])

boxplot(dat.zup)


