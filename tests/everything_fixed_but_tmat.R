library(outbreaker)
library(ape)
options(error=recover)

tm <- matrix(byrow=TRUE,c(0.45,0.45,0.1,0.45,0.45,0.1,0.1,0.1,0.8),ncol=3)
data <- list(n=0)
while(data$n < 10){
data <- simOutbreak(R0=2, infec.curve=c(0,1,1,1),n.hosts=99,duration=10,spatial=FALSE,group.sizes=c(33,33,33),trans.mat=tm)}

collecDates <- data$onset
fake.dna <- data$dna
fake.groups <- data$group
nsize <- data$n
w <- c(rep(1/10,10),rep(0,40))
tre <- data$ances
tre[1] <- -1

res <- outbreaker.parallel(n.runs=4,parallel=TRUE,init.tree=tre,dna=fake.dna, dates=collecDates, idx.dna=c(1:nsize),
                           mut.model=1, spa.model=0,
                           w.dens=w,f.dens=c(0,w),
                           dist.mat=NULL,
                           init.kappa=1, init.mu1=0.01,init.mu2=0.5, init.spa1=NULL,
                           n.iter=6e5, sample.every=500, tune.every=500,
                           burnin=2e5, import.method="genetic",
                           find.import=TRUE,
                           pi.prior1=10, pi.prior2=1, spa1.prior=1,
                           move.mut=FALSE, move.ances=FALSE, move.kappa=FALSE,
                           move.Tinf=FALSE, move.pi=FALSE, move.spa=FALSE, move.Tmat=TRUE,
                           outlier.threshold = 5,
                           quiet=FALSE, res.file.name="chains.txt",
                           tune.file.name="tuning.txt",group.vec=fake.groups)






#res <- outbreaker.parallel(n.runs=4,parallel=TRUE,dna=fake.dna, dates=collecDates, idx.dna=c(1:nsize),
#               mut.model=1, spa.model=0,
#               w.dens=w,f.dens=c(0,w),
#               dist.mat=NULL,
 #              init.tree="seqTrack",
  #             init.kappa=1, init.mu1=0.01,init.mu2=0.5, init.spa1=NULL,
   #            n.iter=6e5, sample.every=500, tune.every=500,
    #           burnin=2e5, import.method="genetic",
     #          find.import=TRUE,
      #         pi.prior1=10, pi.prior2=1, spa1.prior=1,
       #        move.mut=TRUE, move.ances=1, move.kappa=0,
        #       move.Tinf=TRUE, move.pi=TRUE, move.spa=TRUE, move.Tmat=TRUE,
         #      outlier.threshold = 5, max.kappa=rep(1,nsize),
          #     quiet=FALSE, res.file.name="chains.txt",
           #    tune.file.name="tuning.txt",group.vec=fake.groups)

#tre <- get.tTree(res)
#col <- rep("lightgrey",34)
#col[which(data$outbreak$ances != tre$ances)] <- "pink"
#plot(tre, annot="", vertex.color=col)




#par(mfrow=c(2,1))
#plotChains(res)
#plotChains(res,burnin=2e4,type="dens")

par(mfrow=c(3,3))
plot(1:length(res$chains$p_11),res$chains$p_11,type='l')
plot(1:length(res$chains$p_11),res$chains$p_12,type='l')
plot(1:length(res$chains$p_11),res$chains$p_13,type='l')
plot(1:length(res$chains$p_11),res$chains$p_21,type='l')
plot(1:length(res$chains$p_11),res$chains$p_22,type='l')
plot(1:length(res$chains$p_11),res$chains$p_23,type='l')
plot(1:length(res$chains$p_11),res$chains$p_31,type='l')
plot(1:length(res$chains$p_11),res$chains$p_32,type='l')
plot(1:length(res$chains$p_11),res$chains$p_33,type='l')
print(tm)
post.tm <- matrix(byrow=TRUE,ncol=3,c(mean(res$chains$p_11),mean(res$chains$p_12),mean(res$chains$p_13),
                                      mean(res$chains$p_21),mean(res$chains$p_22),mean(res$chains$p_23),
                                      mean(res$chains$p_31),mean(res$chains$p_32),mean(res$chains$p_33)))
print(post.tm)


