library(outbreaker)
library(digest)


# set up parameters #

## GENERATION TIME/SERIAL INTERVAL (VERY SIMILAR, AS SUGGESTED IN PAPER) - UPDATED FIGURES FROM NEW ENG PAPER, MEAN = 13.5, SD = 9.2##
beta <- 13.5/(9.2^2)
alpha <- 13.5*beta
w <- dgamma(seq(0,50,1),shape=alpha,rate=beta)


## COLONIZATION TIME ##
INFEC.CURVE <- w
#INFEC.CURVE <- c(0,1,1,1,1,1)



## TRANSMISSION MATRIX ##
TRANS.MAT <- matrix(ncol=4,byrow=TRUE,
                    c(0.65,0.1,0.15,0.1,
                      0.1,0.6,0.1,0.2,
                      0.05,0.15,0.4,0.4,
                      0.15,0.05,0.4,0.4))

## two big towns (1 & 2)
## two smaller towns (3 & 4)
## 1 closer to 3, 3 and 4 very close, 4 closer to 2, 1 and 2 far apart


## REPRODUCTIVE NUMBERS ##
my.R0 <- 2.1 # VALUE WITHIN CI FOR EACH COUNTRY FROM NEW ENG PAPER


## RATE OF MUTATION ##
transitions <- 2 * (10^-3) *(1/365) # from science mag genomic surveillance paper, gives a substitution rate of ~0.00005 per site per day
transversions <- transitions/2 #acceptable?
dnalength <- 19000 #from http://www.ncbi.nlm.nih.gov/nucleotide/NC_002549


## N.HOSTS / DURATION ##
N.HOSTS <- 200
GRP.SIZES <- c(75,75,25,25)

## MISC ##
SPATIAL <- FALSE
DURATION <- 50


###########################
## SIMULATION LOOP START ##
###########################

for(i in 1:1){

## generate data ##

  #generate data  
  dat <- simOutbreak(R0=my.R0,n.hosts=N.HOSTS,duration=DURATION,spatial=SPATIAL,trans.mat=TRANS.MAT,
                     mu.transi= transitions, mu.transv=transversions,seq.length=dnalength,group.sizes=GRP.SIZES,infec.curve=INFEC.CURVE)
  while(dat$n < 30){
    dat <- simOutbreak(R0=my.R0,n.hosts=N.HOSTS,duration=DURATION,spatial=SPATIAL,trans.mat=TRANS.MAT,
                       mu.transi= transitions, mu.transv=transversions,seq.length=dnalength,group.sizes=GRP.SIZES,infec.curve=INFEC.CURVE)
    print("sim end")
  }






## sort out keys and directories ##
key <- paste(unlist(strsplit(digest(dat),""))[1:10],collapse="")
dir.create(key)
olddir <- getwd()
setwd(key)
cat("Working...", file="ONGOING")


## set up outbreaker params ##
N.ITER = 1e5
BURNIN = 2e4
OUT.DNA <- dat$dna
OUT.GROUPS <- dat$group
OUT.DATES <- dat$onset ## + sample from colon.dens


## PRIORS ##
TMAT.PRIOR.MULT <- 10


              



## run outbreaker ##
print("both")
res <- outbreaker(dna=OUT.DNA,
                  dates = OUT.DATES,
                  mut.model=2,
                  spa.model=0,
                  grp.model=1,
                  w.dens = w,
                  f.dens = w,
                  init.tree = "seqTrack",
                  init.kappa = 1,
                  sample.every=500,tune.every=500,
                  burnin=BURNIN,n.iter=N.ITER,
                  move.Tinf=FALSE,move.kappa=FALSE,
                  group.vec=OUT.GROUPS,
                  tmat.prior.mult=TMAT.PRIOR.MULT)
print("no dna")            
res.nodna <- outbreaker(dates = OUT.DATES,
                        mut.model=0,
                        spa.model=0,
                        grp.model=1,
                        w.dens = w,
                        f.dens = w,
                        init.tree = "seqTrack",
                        init.kappa = 1,
                        sample.every=500,tune.every=500,
                        burnin=BURNIN,n.iter=N.ITER,
                        move.Tinf=FALSE,move.kappa=FALSE,
                        group.vec=OUT.GROUPS,
                        tmat.prior.mult=TMAT.PRIOR.MULT)
print("neither")
res.nodna.nogrps <- outbreaker(dates = OUT.DATES,
                               mut.model=0,
                               spa.model=0,
                               grp.model=0,
                               w.dens = w,
                               f.dens = w,
                               init.tree = "seqTrack",
                               init.kappa = 1,
                               sample.every=500,tune.every=500,
                               burnin=BURNIN,n.iter=N.ITER,
                               move.Tinf=FALSE,move.kappa=FALSE)
print("no groups")
res.nogrps <- outbreaker(dna=OUT.DNA,
                         dates = OUT.DATES,
                         mut.model=2,
                         spa.model=0,
                         grp.model=0,
                         w.dens = w,
                         f.dens = w,
                         init.tree = "seqTrack",
                         init.kappa = 1,
                         sample.every=500,tune.every=500,
                         burnin=BURNIN,n.iter=N.ITER,
                         move.Tinf=FALSE,move.kappa=FALSE)


## shape output ##

file.remove("ONGOING")
save(res,res.nodna,res.nodna.nogrps,res.nogrps,dat,key, file=paste(key,"RData",sep="."))
setwd(olddir)


}

#########################
## SIMULATION LOOP END ##
#########################

