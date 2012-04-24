
##################################
## function to get all ancestries
##################################
getAnces <- function(x){
    ## get ancestries for one infection
    f1 <- function(x, id){
        ## check for auto-infection
        if(any(x == 1:length(x))) {
            cat("\nAuto-infection detected:\n")
            print(which(x == 1:length(x)))
            stop()
        }

        ances <- id
        while(ances[1]>0){
            ances <- c(x[ances[1]],ances)
            ## check for cycles
            if(ances[1] %in% ances[-1]){
                cat("\nCycles detected.\n")
                break
            }
        }
        return(ances)
    }

    res <- lapply(1:length(x), function(i) f1(x, i))
    return(res)
}



################################
## function to check ancestries
################################
checkAnces <- function(listAnces){
    ## returns TRUE if OK
    isOK <- function(ances){
        if(any(table(ances)>1)) return(FALSE)
        if(ances[1]==-100 & length(ances)>2) return(FALSE)
        return(TRUE)
    }

    res <- sapply(listAnces, isOK)
    return(res) 
}



##############
## application
##############
x <- scan("IndexInfector.txt")
x[x>0] <- x[x>0]+1 # fix for C->R vector indexing

toto <- getAnces(x)
toto # all ancestries
checkAnces(toto)
toto[!checkAnces(toto)] # problematic ones
