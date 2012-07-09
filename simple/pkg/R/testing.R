
test.outbreaker <- function(){
    if(!require(adegenet)) stop("adegenet is required")

    w <- c(0.1,1,4,3,2,1)
    w <- w/sum(w)
    genTime <- function() {
        return(sample(0:(length(w)-1), 1, prob=w))
    }

    dat <- haploGen(seq.length=1e4, mu=1e-4, t.max=20,
                   gen.time=genTime,
                   repro=function(){round(rnorm(1,2,1))}, max.nb.haplo=1e5,
                   geo.sim=FALSE)

    outbreaker()

}
