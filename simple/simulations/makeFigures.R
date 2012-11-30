
## GET POOLED RESULTS ##
library(adegenet)

source("/home/thibaut/dev/outbreaker/outbreaker-code/simple/simulations/poolResults.R")
x <- poolResults()


## get ordered levels for types of simulation
temp <- sort(unique(gsub("-fixedMu","",x$type)))
lev <- rep(temp, each=2)
change.idx <- seq(2,by=2,length=length(temp))
lev[change.idx] <- paste(lev[change.idx],"fixedMu",sep="-")

type <- factor(x$type, levels=lev) # new factor
all(as.character(type)==x$type) # sanity check - must be TRUE


###############################
## FIGURE 1: ANCESTRY RESULTS
###############################
temp <- sort(unique(gsub("-fixedMu","",x$type)))
type.col <- rep(fac2col(temp, col.pal=rainbow),each=2)

## pdf ##
pdf("figures/figure1.pdf")
par(mar=c(6.1,4.1,1.1,1.1))

## boxplot
boxplot(x$prop.ances.ok~type, col=type.col,xaxt="n", ylab="Proportion of correct ancestries in consensus tree", cex.lab=1.2)

## lines
abline(v=seq(0.5,by=2,le=20),col="lightgrey")

## axes
ax.at <- seq(1.5, by=2, length=length(temp))
axis(side=1, at=ax.at, label=temp,las=3)
dev.off()


## svg ##
svg("figures/figure1.svg")
par(mar=c(6.1,4.1,1.1,1.1))

## boxplot
boxplot(x$prop.ances.ok~type, col=type.col,xaxt="n", ylab="Proportion of correct ancestries in consensus tree", cex.lab=1.2)

## lines
abline(v=seq(0.5,by=2,le=20),col="lightgrey")

## axes
ax.at <- seq(1.5, by=2, length=length(temp))
axis(side=1, at=ax.at, label=temp,las=3)
dev.off()









#########################################
## FIGURE S1: MEAN INFECTION DATE ERRORS
#########################################
temp <- sort(unique(gsub("-fixedMu","",x$type)))
type.col <- rep(fac2col(temp, col.pal=rainbow),each=2)

## pdf ##
pdf("figures/figureS1.pdf")
par(mar=c(6.1,4.1,1.1,1.1))

## boxplot
boxplot(x$merr.date~type, col=type.col,xaxt="n", ylab="Mean error in estimated infection date", cex.lab=1.2)

## lines
abline(v=seq(0.5,by=2,le=20),col="lightgrey")

## axes
ax.at <- seq(1.5, by=2, length=length(temp))
axis(side=1, at=ax.at, label=temp,las=3)
dev.off()

## svg ##
svg("figures/figureS1.svg")
par(mar=c(6.1,4.1,1.1,1.1))

## boxplot
boxplot(x$merr.date~type, col=type.col,xaxt="n", ylab="Mean error in estimated infection date", cex.lab=1.2)

## lines
abline(v=seq(0.5,by=2,le=20),col="lightgrey")

## axes
ax.at <- seq(1.5, by=2, length=length(temp))
axis(side=1, at=ax.at, label=temp,las=3)
dev.off()







#####################################
## FIGURE S2: IMPORTED CASE RESULTS
#####################################
temp <- sort(unique(gsub("-fixedMu","",x$type)))
type.col <- rep(fac2col(temp, col.pal=rainbow),each=2)

## pdf ##
pdf("figures/figureS2.pdf")
par(mar=c(6.1,4.1,1.1,1.1))

## boxplot
boxplot(x$prop.imp.ok~type, col=type.col,xaxt="n", ylab="Proportion of imported cases detected", cex.lab=1.2)

## lines
abline(v=seq(0.5,by=2,le=20),col="lightgrey")

## axes
ax.at <- seq(1.5, by=2, length=length(temp))
axis(side=1, at=ax.at, label=temp,las=3)
dev.off()

## svg ##
svg("figures/figureS2.svg")
par(mar=c(6.1,4.1,1.1,1.1))

## boxplot
boxplot(x$prop.imp.ok~type, col=type.col,xaxt="n", ylab="Proportion of imported cases detected", cex.lab=1.2)

## lines
abline(v=seq(0.5,by=2,le=20),col="lightgrey")

## axes
ax.at <- seq(1.5, by=2, length=length(temp))
axis(side=1, at=ax.at, label=temp,las=3)
dev.off()






##########################################
## FIGURES FOR HEREROGENEITY SIMULATIONS
##########################################
source("analyseHetero.R")

dat.2g <- cbind.data.frame(trueR[1:50,], dat.2g[1:50,c(2,3,5,6)])
dat.zup <- dat.zup[1:50,]



pdf("figures/figure2.pdf")

## margins
par(mar=c(5.1,4.1,2.1,4.1))

## boxplot for 2-group simulation
lab.2g <- rep(c("Low R", "High R"),3)
col.2g <- rep(c("grey","deepskyblue3","red3"),each=2)
boxplot(dat.2g, names=lab.2g, col=col.2g, ylab="Effective reproduction number (R)",cex.lab=1.2, xlim=c(0,10),las=3)

## separation
abline(v=7, lty=2, lwd=2)

## boxplot for super-spreaders
ousr <- par("usr")
par(usr=c(ousr[1:2],-0.05,1.05))

boxplot(dat.zup,at=8:9,add=TRUE, col=c("deepskyblue3","red3"), yaxt="n",xaxt="n")

## add axis
axis(side=4,srt=180)
mtext(side=4, "Proportion of super-spreader detected", cex=1.2, line=3,srt=180)

## legend
par(xpd=TRUE)
leg.xy <- list(x=6.5, y=par("usr")[3])
legend(leg.xy, fill=c("grey","deepskyblue3","red3"), legend=c("actual values", "model with DNA","model without DNA"),bg=transp("white",.7),cex=1)

dev.off()



svg("figures/figure2.svg")

## margins
par(mar=c(5.1,4.1,2.1,4.1))

## boxplot for 2-group simulation
lab.2g <- rep(c("Low R", "High R"),3)
col.2g <- rep(c("grey","deepskyblue3","red3"),each=2)
boxplot(dat.2g, names=lab.2g, col=col.2g, ylab="Effective reproduction number (R)",cex.lab=1.2, xlim=c(0,10),las=3)

## separation
abline(v=7, lty=2, lwd=2)

## boxplot for super-spreaders
ousr <- par("usr")
par(usr=c(ousr[1:2],-0.05,1.05))

boxplot(dat.zup,at=8:9,add=TRUE, col=c("deepskyblue3","red3"), yaxt="n",xaxt="n")

## add axis
axis(side=4,srt=180)
mtext(side=4, "Proportion of super-spreader detected", cex=1.2, line=3,srt=180)

## legend
par(xpd=TRUE)
leg.xy <- list(x=6.5, y=par("usr")[3])
legend(leg.xy, fill=c("grey","deepskyblue3","red3"), legend=c("actual values", "model with DNA","model without DNA"),bg=transp("white",.7),cex=1)

dev.off()











## ## general results - consensus ancestry OK
## qplot(prop.ances.ok, data=x, color=type, geom="density")

## ## general results - actual ancestor support
## qplot(msup.ances, data=x, color=type, geom="density")

## ## support for actual kappa
## qplot(msup.kappa, data=x, color=type, geom="density")

## ## support for mu1
## qplot(mu1.ok, data=x, fill=type, geom="bar")

## ## support for mu2
## qplot(mu2.ok, data=x, fill=type, geom="bar")

## ## error about the date
## qplot(merr.date, data=x, color=type, geom="density")

## ## support for mu2
## qplot(prop.imp.ok, data=x, fill=type, geom="bar")

## ## imported cases
## qplot(prop.imp.ok, data=x, fill=type, geom="bar")

## ## time vs n
## qplot(n, time/60, data=x, color=type, geom=c("point"), ylab="time (minutes)")
