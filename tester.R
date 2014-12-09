library(igraph)
library(outbreaker)	
library(adegenet)

dm <- matrix(byrow=TRUE,ncol=2,c(0,1,1,0))
res <- simOutbreak(R0=c(2,2.5,1),infec.curve=c(0,1,1,1),n.hosts=40,spatial=FALSE,group.sizes=c(20,20),trans.mat=dm,duration=5)

plot(as.igraph.simOutbreak(res,vertex.col="group"))
legend(x="topleft",col=num2col(unique(res$group),col.pal=funky),legend=c(1:2),pch=16)
