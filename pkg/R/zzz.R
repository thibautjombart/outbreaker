.onAttach <- function(libname, pkgname){
    pkg.version <- packageDescription("outbreaker", fields = "Version")

    startup.txt <- paste("\n   === outbreaker", pkg.version, "is loaded ===   \n\n")
    startup.txt <- paste(startup.txt, "Questions/documentation -> check the R-epi project: \nhttp://sites.google.com/site/therepiproject", sep="\n")
    packageStartupMessage(startup.txt)
}
