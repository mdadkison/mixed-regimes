niters <-1000
big <- array(dim=niters)
outname <- "mixedResults/BigSmsy.txt"          
write("Fraction of Smsy estimates greater than one million",file=outname)


for (man in c("Precise","Imprecise")) {
  for(harv in c("Overh","Normh","Underh")) {
    for(sh in c("Same","Down")) {
      #browser()
      if(sh=="Same") {
        for(st in c("Histart","Lostart")) {
          star <- st
          fname <- paste("mixedResults/R",man,harv,sh,star,"Sims.txt")
          fname <- gsub(" ", "", fname, fixed = TRUE)
          xfile <- read.table(fname,header=TRUE)
          big <- ifelse(xfile$Smsy>200000,1,0)
          frac_big <- sum(big/niters)
          cat(fname,frac_big,"\n")
          write(c(fname,frac_big),file=outname,append=TRUE)
          
        } 
      }else {
          star <-"Histart"
          fname <- paste("mixedResults/R",man,harv,sh,star,"Sims.txt")
          fname <- gsub(" ", "", fname, fixed = TRUE)
          xfile <- read.table(fname,header=TRUE)
          big <- ifelse(xfile$Smsy > 1000000,1,0)
          frac_big <- sum(big/niters)
          cat(fname,frac_big,"\n")
          write(c(fname,frac_big),file=outname,append=TRUE)
          
      }
        
        
    } # sh
  } # harv
} # man

