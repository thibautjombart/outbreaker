.onAttach <- function(libname, pkgname){
    pkg.version <- packageDescription("outbreaker", fields = "Version")

    startup.txt <- paste("\n   === outbreaker", pkg.version, "is loaded ===   \n\n")
    packageStartupMessage(startup.txt)
}
