.First.lib <- function (lib, pkg){
    library.dynam("outbreaker", pkg, lib)
    pkg.version <- packageDescription("outbreaker", fields = "Version")

    startup.txt <- paste("   ==========================\n    outbreaker", pkg.version, "is loaded\n   ==========================\n\n")

    packageStartupMessage(startup.txt)
}
