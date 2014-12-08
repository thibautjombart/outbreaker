library(igraph)
library(outbreaker)	

dm <- matrix(byrow=TRUE,ncol=3,c(0.998,0.001,0.001,0.001,0.998,0.001,0.001,0.001,0.998))
res <- simOutbreak(R0=c(2,2.5,1),infec.curve=c(0,1,1,1),n.hosts=150,spatial=FALSE,group.sizes=c(50,50,50),trans.mat=dm,duration=8)

plot(as.igraph.simOutbreak(res,vertex.col="group"))
legend(x="topleft",col=num2col(unique(res$group),col.pal=heat.colors),legend=c(1:3),pch=1)
