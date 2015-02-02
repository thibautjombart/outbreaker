library(outbreaker)
library(ape)
library(digest)
library(adegenet)
library(EpiEstim)
options(error=recover)


# set transmission matrix
tm <- matrix(byrow=TRUE,c(1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3),ncol=3)


for(i in 1:5){
  
  #generate data  
  dat <- simOutbreak(R0=2, infec.curve=c(0,1,1,1),n.hosts=99,duration=10,spatial=FALSE,group.sizes=c(33,33,33),trans.mat=tm)
  while(dat$n < 30){
    dat <- simOutbreak(R0=2, infec.curve=c(0,1,1,1),n.hosts=99,duration=10,spatial=FALSE,group.sizes=c(33,33,33),trans.mat=tm)
  }
  
  #create key, sort out directories
  key <- paste(unlist(strsplit(digest(dat),""))[1:10],collapse="")
  dir.create(key)
  olddir <- getwd()
  setwd(key)
  cat("Working...", file="ONGOING")
  
  
  
  #set up outbreaker
  collecDates <- dat$onset
  fake.dna <- dat$dna
  fake.groups <- dat$group
  nsize <- dat$n
  
  # create gentime
  max.w <- max(dat$onset)-min(dat$onset)
  w <- c(seq(from=0,to=max.w/2),seq(from=max.w/2,to=0))
  w <- w + dnorm(w,mean=max.w/2,sd=2)
  w <- w/sum(w)
  
  BURNIN = 2e4
  
  
  print(paste("run number ", i))
  res <- outbreaker(dna=fake.dna, dates=collecDates,
                    mut.model=2, spa.model=0,
                    w.dens=w,
                    init.kappa=1,
                    n.iter=1e5, sample.every=500, tune.every=500,
                    burnin=BURNIN,
                    find.import=TRUE,
                    move.kappa=FALSE,
                    move.Tinf=FALSE, 
                    res.file.name="chains.txt",
                    tune.file.name="tuning.txt",group.vec=fake.groups)
  
  #shape output
  tre <- get.tTree(res,burnin=BURNIN)
  chains <- res$chains[res$chains$step>BURNIN,,drop=FALSE]
  
  
  ubmat <- matrix(ncol=3,byrow=TRUE,
                  c(quantile(res$chains$p_11,0.975),
                    quantile(res$chains$p_12,0.975),
                    quantile(res$chains$p_13,0.975),
                    quantile(res$chains$p_21,0.975),
                    quantile(res$chains$p_22,0.975),
                    quantile(res$chains$p_23,0.975),
                    quantile(res$chains$p_31,0.975),
                    quantile(res$chains$p_32,0.975),
                    quantile(res$chains$p_33,0.975))
  )
  
  
  lbmat <- matrix(ncol=3,byrow=TRUE,
                  c(quantile(res$chains$p_11,0.025),
                    quantile(res$chains$p_12,0.025),
                    quantile(res$chains$p_13,0.025),
                    quantile(res$chains$p_21,0.025),
                    quantile(res$chains$p_22,0.025),
                    quantile(res$chains$p_23,0.025),
                    quantile(res$chains$p_31,0.025),
                    quantile(res$chains$p_32,0.025),
                    quantile(res$chains$p_33,0.025))
  )
  
  
  
  #record success/failure
  #deal with files
  file.remove("ONGOING")
  save(res,dat,key,tre,chains,ubmat,lbmat, file=paste(key,"RData",sep="."))
  setwd(olddir)
  
  
  
}