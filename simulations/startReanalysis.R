source("~/dev/outbreaker/outbreaker-code/simple/simulations/poolResults.R")

source("~/dev/outbreaker/outbreaker-code/simple/simulations/reanalyseSimul.1.0.R")


source("../reanalyseSimul.1.0.R")
sapply(dir()[2:10],reanalyseSimul)

source("../reanalyseSimul.1.0.R")
sapply(dir()[21:30],reanalyseSimul)

source("../reanalyseSimul.1.0.R")
sapply(dir()[31:40],reanalyseSimul)

source("../reanalyseSimul.1.0.R")
sapply(dir()[41:50],reanalyseSimul)
