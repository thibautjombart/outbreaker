
in.files <- list.files(path=getwd())
olddir <- getwd()
in.files <- in.files[-grep(in.files,pattern=".R")]
winmat <- matrix(ncol=3,rep(0,9))
rowSkip <- c(2,3,1)

for(current in in.files){
	setwd(paste(current))
	fname <- paste(current,"RData",sep=".")
	load(fname)

  for(i in 1:3){
    for(j in 1:3){
      if(lbmat[i,j]<1/3 & ubmat[i,j] > 1/3){winmat[i,j] <- winmat[i,j]+ 1}
    }
  }
  
  rm(fname)
	setwd(olddir)
}
print("Proportion of parameter credible intervals that include the true value")
print(winmat/length(in.files))