
library(outbreaker)
set.seed(2)

## TEST USING TOYOUTBREAK ##
data(toyOutbreak)
attach(toyOutbreak)
plot(dat)

##res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4)
res <-  outbreaker.parallel(n.runs=6, dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4)
plotChains(res)

temp1 <- dat$ances
temp2 <- get.tTree(res)$ances
temp1[is.na(temp1)] <- 0
temp2[is.na(temp2)] <- 0
df <- data.frame(true=temp1, recons=temp2, diff=ifelse(temp1!=temp2, "!", " "))
df
mean(temp1==temp2)

areDifferent <- which(temp1 != temp2)
discrep <- lapply(areDifferent, function(i) table(res$chains[[paste("alpha",i,sep="_")]], res$chains[[paste("kappa",i,sep="_")]]))
names(discrep) <- paste("ancestry of", areDifferent, ":", apply(df[areDifferent,], 1, paste, collapse=" "))
discrep



w <- c(0, 0.5, 1, 0.75)

## NON-SPATIAL SIMULATION ##
## this may generate an error if outbreak doesn't take off
dat <- simOutbreak(R0 = 2, infec.curve = w, n.hosts = 100, spatial=FALSE, mu.transi=.5e-4)[1:30]
##collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)
collecDates <- dat$onset + 2
plot(dat)

## test parallel
#res <-  outbreaker.parallel(n.runs=2, dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5)

## run outbreaker
##res <-  outbreaker(dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4, spa.model=0)
res <-  outbreaker.parallel(3,dna=dat$dna, dates=collecDates,w.dens=w, n.iter=5e4, spa.model=0)

plotChains(res)

plot(dat)
temp1 <- dat$ances
temp2 <- get.tTree(res)$ances
temp1[is.na(temp1)] <- 0
temp2[is.na(temp2)] <- 0
df <- data.frame(true=temp1, recons=temp2, diff=ifelse(temp1!=temp2, "!", " "))
df
mean(temp1==temp2)

areDifferent <- which(temp1 != temp2)
discrep <- lapply(areDifferent, function(i) table(res$chains[[paste("alpha",i,sep="_")]], res$chains[[paste("kappa",i,sep="_")]]))
names(discrep) <- paste("ancestry of", areDifferent, ":", apply(df[areDifferent,], 1, paste, collapse=" "))
discrep


## SPATIAL SIMULATION ##
## test spatial model 1 (exponential)
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)

## this may generate an error if outbreak doesn't take off
dat <- simOutbreak(R0 = 1.5, infec.curve = w, n.hosts = 100, spatial=TRUE)[1:35]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)
D <- as.matrix(dist(dat$xy))

res <-  outbreaker.parallel(3,dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=5e4, spa.model=1, find.import=FALSE)
plotChains(res)

temp1 <- dat$ances
temp2 <- get.tTree(res)$ances
temp1[is.na(temp1)] <- 0
temp2[is.na(temp2)] <- 0
df <- data.frame(true=temp1, recons=temp2, diff=ifelse(temp1!=temp2, "!", " "))
df
mean(temp1==temp2)

areDifferent <- which(temp1 != temp2)
discrep <- lapply(areDifferent, function(i) table(res$chains[[paste("alpha",i,sep="_")]], res$chains[[paste("kappa",i,sep="_")]]))
names(discrep) <- paste("ancestry of", areDifferent, ":", apply(df[areDifferent,], 1, paste, collapse=" "))
discrep


dist.inf <- na.omit(sapply(1:dat$n, function(i) D[dat$id[i], dat$ances[i]]))
mean(dist.inf)

par(mfrow=c(2,2))
plotChains(res)
plotChains(res, what="spa1")
plotChains(res, what="spa1", type="dens")
abline(v=mean(dist.inf))
hist(dist.inf , col="grey",
     main="Distance of infectious contacts", xlab="Distance infector-infected")




## bigger example ##
library(outbreaker)
w <- c(0, 0.5, 1, 0.75)

## this may generate an error if outbreak doesn't take off
dat <- simOutbreak(R0 = 1.5, infec.curve = w, n.hosts = 200, spatial=TRUE)[1:50]
collecDates <- dat$onset + sample(0:3, size=length(dat$onset), replace=TRUE, prob=w)
D <- as.matrix(dist(dat$xy))

res <-  outbreaker.parallel(3, dna=dat$dna, dates=collecDates,w.dens=w, dist.mat=D, n.iter=1e5, spa.model=1)

temp1 <- dat$ances
temp2 <- get.tTree(res)$ances
temp1[is.na(temp1)] <- 0
temp2[is.na(temp2)] <- 0
df <- data.frame(true=temp1, recons=temp2, diff=ifelse(temp1!=temp2, "!", " "))
df
mean(temp1==temp2)

areDifferent <- which(temp1 != temp2)
discrep <- lapply(areDifferent, function(i) table(res$chains[[paste("alpha",i,sep="_")]], res$chains[[paste("kappa",i,sep="_")]]))
names(discrep) <- paste("ancestry of", areDifferent, ":", apply(df[areDifferent,], 1, paste, collapse=" "))
discrep

dist.inf <- na.omit(sapply(1:dat$n, function(i) D[dat$id[i], dat$ances[i]]))
mean(dist.inf)

par(mfrow=c(2,2))
plotChains(res)
plotChains(res, what="spa1")
plotChains(res, what="spa1", type="dens")
abline(v=mean(dist.inf))
hist(dist.inf , col="grey",
     main="Distance of infectious contacts", xlab="Distance infector-infected")




### test missing sequences ###
library(outbreaker)
set.seed(2)

## TEST USING TOYOUTBREAK ##
data(toyOutbreak)
attach(toyOutbreak)
dat <- dat[1:15]
plot(dat)

seq.for.cases <- c(1,2,3,8,15,4,9,7,11)

res <-  outbreaker(dna=dat$dna[seq.for.cases,], idx.dna=seq.for.cases, dates=collecDates,w.dens=w, n.iter=5e4)
