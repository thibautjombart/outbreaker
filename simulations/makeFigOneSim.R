###########################################################
## SCRIPT USED TO GENERATE FIGURES FOR A GIVEN SIMULATION
###########################################################

## LOAD PACKAGES ##
library(outbreaker)
library(adegenet)
library(ape)


## CLEAN WORKSPACE ##
rm(list=ls())


## FIND DATA AND RESULTS ##
sim <- "4d722dc202"

## find files locations
files <- grep(sim, dir(recursive=TRUE), value=TRUE)
f.in <- grep("in.csv", files, value=TRUE)
f.out <- grep("out.csv", files, value=TRUE)
f.RDat <- grep("RData", files, value=TRUE)

## read files, load objects
inputs <- read.csv(f.in)
outputs <- read.csv(f.out)
load(f.RDat)
ls()




## MAKE FIGURES ##
BURNIN <- 2e4

pdf("figures/oneSim-post.pdf")
par(mfrow=c(2,1), mar=c(4,3,3,2))
plot.chains(res, col=funky(6), main="MCMC - posterior chains")
box("fig")
plot.chains(res, col=funky(6), type="dens", main="MCMC - posterior density (without burnin)", omit=BURNIN)
box("fig")
dev.off()


## get consensus tree
x <- get.TTree.simple(res, burnin=BURNIN)
mean(x$ances==dat$ances,na.rm=TRUE)


## ## INFECTIOUSNESS CURVES, FIRST 6 CASES ##
## pdf("figures/oneSim-epicurves.pdf")
## epicurves(x[1:6,], coef=2, cex.lab=1.25)
## dev.off()



## DATA VS RECONSTRUCTION ##
pdf("figures/oneSim-data-recons.pdf")

## color for the vertices
v.col <- rep("lightblue",length(x$ances))
notOk <- which(x$ances!=dat$ances)
if(length(notOk)>0) v.col[notOk] <- "red"
par(mfrow=c(1,2), mar=c(.5,.5,4,.5))
plot(dat,main="Simulated data", edge.color="black", annot="")

## reconstructed graph
g <- as.igraph(x)
arr.col <- rep("lightgrey", length(E(g)$prob))
arr.col[E(g)$prob>0.5] <- "black"
plot(g ,main="Reconstructed outbreak", vertex.color=v.col,
     edge.color=arr.col, layout=attr(g, "layout"), edge.label="", edge.width=2)
mtext(side=3, text="(erroneous inferences in red)")
legend("bottomleft", lty=1, lwd=4, col=c("grey", "black"), leg=c("support < 50%", "support >= 50%"))

dev.off()





## FIGURE FOR TINF ##
pdf("figures/oneSim-Tinf.pdf")

toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
colnames(Tinf) <- sub("Tinf_", "case ",colnames(Tinf))
boxplot(Tinf, col=funky(20), horizontal=TRUE, las=1, xlab="Time", main="Inference of infection dates")
mtext(side=3, text="(X = actual date)")
points(dat$dates,1:length(dat$dates), col="black", pch="x",cex=2.5)
points(dat$dates,1:length(dat$dates), col=funky(20), pch="x",cex=2)

dev.off()







## ## check results ##

## ## check ancestries
## x <- get.TTree.simple(res, burn=BURNIN)
## temp <- x
## temp$ances[is.na(temp$ances)] <- 0
## temp2 <- dat$ances
## temp2[is.na(temp2)] <- 0
## temp2[1] <- NA
## mean(temp$ances==temp2,na.rm=TRUE)

## v.col <- rep("lightblue",length(x$ances))
## notOk <- which(temp$ances!=temp2)
## if(length(notOk)>0) v.col[notOk] <- "red"
## par(mfrow=c(1,1))
## plot(dat,main="data", vertex.color=v.col)
## x11();
## plot(x,main="reconstruction (red=wrong ancestry)", vertex.color=v.col)


## ## check frequency of external infections
## alpha <- res$chains[,grep("alpha", names(res$chains))]
## barplot(apply(alpha,2, function(e) mean(e==0)), main="Proportion of inferred external case")


## ## check mutation rates
## par(mfrow=c(2,1))
## plot.chains(res, "mu1",type="dens", omit=2e4)
## abline(v=1e-4, col="blue")
## plot.chains(res, "mu2",type="dens", omit=2e4)
## abline(v=1e-4, col="blue")



## ## check kappa
## v.col <- rep("lightblue",length(x$ances))
## notOk <- which(x$n.gen!=1)
## if(length(notOk)>0) v.col[notOk] <- "red"
## plot(dat,main="data", vertex.color=v.col)
## x11();
## plot(x,main="reconstruction (red=wrong kappa)", vertex.color=v.col, annot="n.gen")



## ## check Tinf
## toKeep <- grep("Tinf",names(res$chains))
## Tinf <- res$chains[res$chains$step>BURNIN, toKeep]
## boxplot(Tinf, col="grey")
## points(dat$dates, col="red", pch="x",cex=2)



## ## check Pi
## plot.chains(res, "pi", omit=BURNIN, type="de")
## abline(v=1)



############################################





