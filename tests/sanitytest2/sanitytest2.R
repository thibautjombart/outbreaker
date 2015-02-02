library(outbreaker)
library(ape)
library(digest)
library(adegenet)
library(EpiEstim)
options(error=recover)

win=0
# set transmission matrix
tm <- matrix(byrow=TRUE,c(0,1,0,0,0,1,1,0,0),ncol=3)


for(i in 1:2){

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
within1 = TRUE
collecDates <- dat$onset
fake.dna <- dat$dna
fake.groups <- dat$group
nsize <- dat$n
w <- c(0,0.25,0.75,1,0.5)
BURNIN = 2e4


print(paste("run number ", i))
res <- outbreaker(dna=fake.dna, dates=collecDates, idx.dna=c(1:nsize),
                           mut.model=2, spa.model=0,
                           w.dens=w,f.dens=c(0,w),
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






#check output
    if(mean(res$chains$p_11)-1.96*sd(res$chains$p_11) > 0){within1==FALSE}
    if(mean(res$chains$p_12)+1.96*sd(res$chains$p_12) < 1){within1==FALSE}
    if(mean(res$chains$p_13)-1.96*sd(res$chains$p_13) > 0){within1==FALSE}
    if(mean(res$chains$p_21)-1.96*sd(res$chains$p_21) > 0){within1==FALSE}
    if(mean(res$chains$p_22)-1.96*sd(res$chains$p_22) > 0){within1==FALSE}
    if(mean(res$chains$p_23)+1.96*sd(res$chains$p_23) < 1){within1==FALSE}
    if(mean(res$chains$p_31)+1.96*sd(res$chains$p_31) < 1){within1==FALSE}
    if(mean(res$chains$p_32)-1.96*sd(res$chains$p_32) > 0){within1==FALSE}
    if(mean(res$chains$p_33)-1.96*sd(res$chains$p_33) > 0){within1==FALSE}


#record success/failure
if(within1==TRUE){win = win +1}


#deal with files
file.remove("ONGOING")
save(res,dat,key,tre,chains,win, file=paste(key,"RData",sep="."))
setwd(olddir)



}