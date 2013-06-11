## TRY TO SEE HOW THE METHOD BEHAVES ##

############################################
## DATA SIMULATION ##

## load packages and data
library(outbreaker)
library(adegenet)
library(ape)


load("Robjects/data4.RData")
plot(dat, main="Data")
############################################




############################################
## ESTIMATE EVERYTHING  ##
## run outbreaker
## system.time(res <- outbreaker.parallel(n.runs=6, dna=dat$dna, dates=collecDates, w.dens=w, init.tree="seqTrack", n.iter=2e5))
## save(res, file="Robjects/talkCF-res.RData")

load("Robjects/talkCF-res.RData")



## PLOT POSTERIOR ##
pdf("figures/talkCF-post.pdf")
par(mfrow=c(2,1), mar=c(4,3,3,2))
plot.chains(res, col=funky(6), main="MCMC - posterior chains")
box("fig")
plot.chains(res, col=funky(6), type="dens", main="MCMC - posterior density")
box("fig")
dev.off()


## get consensus tree
x <- get.TTree.simple(res)
mean(x$ances==dat$ances,na.rm=TRUE)



## INFECTIOUSNESS CURVES, FIRST 6 CASES ##
pdf("figures/talkCF-epicurves.pdf")
epicurves(x[1:6,], coef=2, cex.lab=1.25)
dev.off()



## DATA VS RECONSTRUCTION ##
pdf("figures/talkCF-data-recons.pdf")

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
pdf("figures/talkCF-Tinf.pdf")

toKeep <- grep("Tinf",names(res$chains))
Tinf <- res$chains[res$chains$step>1e5, toKeep]
colnames(Tinf) <- sub("Tinf_", "case ",colnames(Tinf))
boxplot(Tinf, col=funky(20), horizontal=TRUE, las=1, xlab="Time", main="Inference of infection dates")
mtext(side=3, text="(X = actual date)")
points(dat$dates,1:20, col="black", pch="x",cex=2.5)
points(dat$dates,1:20, col=funky(20), pch="x",cex=2)

dev.off()





