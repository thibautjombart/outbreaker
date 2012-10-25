

## CHECK TRANSITIONS ##
## only transitions, 1 / generations on average ##
dat <- simOutbreak(R0=5, rate.import.case=0, infec.curve=c(0,1), n.hosts=100,mu.transi=1e-4, mu.transv=0)

plot(dat, annot="dist")


## check distribution of mutations
## empirical
f <- table(dat$nmut)/length(dat$nmut)
plot(f, col="blue")

## expected
points(0:10, dpois(0:10, 1), col="red")

## check the nature of the mutations
temp <- as.character(dat$dna[,(seg.sites(dat$dna))])
apply(temp,2,unique)
transi <- table(apply(temp,2,function(v) paste(sort(unique(v)),collapse="-")))
transi







## CHECK TRANSVERSIONS ##
## only transitions, 1 / generations on average ##
dat <- simOutbreak(R0=5, rate.import.case=0, infec.curve=c(0,1), n.hosts=100,mu.transi=0, mu.transv=1e-4)

plot(dat, annot="dist")


## check distribution of mutations
## empirical
f <- table(dat$nmut)/length(dat$nmut)
plot(f, col="blue")

## expected
points(0:10, dpois(0:10, 1), col="red")

## check the nature of the mutations
temp <- as.character(dat$dna[,(seg.sites(dat$dna))])
apply(temp,2,unique)
transv <- table(apply(temp,2,function(v) paste(sort(unique(v)),collapse="-")))
transv
