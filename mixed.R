#Simulates stock-recruitment data and escapement goal estimation when there are
#   autocorrelated environmental effects and MLE estimation
#Question: Is there a tendency for escapement goal estimates to track environmental fluctuations?
#Question: Does the management response result in a downward trend in goals?


#setwd("C:/Users/adkison/Dropbox/My Documents/Active Research/SR review/Habitat") #desktop
#setwd("C:/Users/Milo/Dropbox/My Documents/Active Research/SR review/Habitat") #laptop
library("ggplot2")
library(stats)
source("ShiftyFun.R")

#dummy loop
for (idummy in 1:1) {
  niters <- 1000
  
  ############## simulate data ############################################
  # stock-recruitment parameters for R = aS/(1+aS/b) with log-normal error
  #SRmodel <- "B" #Beverton-Holt
  SRmodel <- "R" #Ricker (changed model to aSexp(1-S/b))
  nyears <- 20
  maxRpS <- exp(1.5)
  capacity <- 100000.
  #maxRpS <- exp(1.5)*0.5
  #capacity <- 100000*0.5
  RCV <- 0.5
  SCV <- 0.5
  phi <- 0.0 #Connors analysis gives 0.43 median for sockeye stocks
  if(SRmodel=="B") {
    alpha <- maxRpS
    beta <- capacity*maxRpS/(maxRpS-1)
    trueSmsy <- BSmsy(alpha,beta)
  } else {
    alpha <- log(maxRpS)
    beta <- capacity
    trueSmsy <- RSmsyExact(c(alpha,beta))$Smsy
  }
  S <- trueSmsy
  
  ########################## MANUAL CHANGES TO SCENARIO IN THIS SECTION ###################################
  leapyear <- 11
  #leapa <- log(2) #shift up
  #leapb <- 2
  #leapa <- log(0.5) #shift down
  #leapb <- 0.5
  leapa <- log(1.0) #no shift
  leapb <- 1.0
  leapduration <- 10
  if(SRmodel=="B") {
    leapSmsy <- BSmsy(alpha*leapa,beta*leapb)
  } else {
    leapSmsy <- RSmsyExact(c(alpha+leapa,beta*leapb))$Smsy
  }
  
  # management options
  #management <- "Precise"
  management <- "Imprecise"
  
  if(beta==100000) {
    hilo <- "Histart"
  } else {
    hilo <- "Lostart"
  }
    
  
  if(management=="precise") {
    Emin <- trueSmsy #constant escapement goal policy
    hrate <- 1
  } else {
    Emin <- 0.8*trueSmsy #constant escapement goal policy
    hrate <- 0.8
  }
  
  #harvHist <- "Normh"
  #harvHist <- "Overh"
  harvHist <- "Underh"
  
  if(harvHist=="Overh") {
    Emin=Emin/2
  } else if (harvHist=="Underh") {
    Emin=Emin*1.5
  }
  
  CVE <- 0.2
  
  if (SCV>0) {
    serr <- "SpErr"
  } else {
    serr <- ""
  }
  
  if(leapb>1){
    leap <- "Up"} else {
      if(leapb<1) {
        leap="Down"
      } else {
        leap="Same"
      }
    }
#  outname <- paste("mixedResults/",SRmodel,management,harvHist,leap,hilo)
  outname <- paste("mixedSError/",SRmodel,management,harvHist,leap,hilo,serr)
  #outname = "junk"
  outname <- gsub(" ", "", outname, fixed = TRUE)
  ############################ end of manual changes ###############################################################
  
  
  ############
  NewGoal <- rep(FALSE,nyears)
  NewGoal[c(20)] <- TRUE
  #NewGoal[c(12,15,18,21,24,27,30,39,49)] <- TRUE
  max_goal_change <- 0.5 #maximum allowed change in escapement goal
  
  spawners <- array(dim=nyears+1)
  spawners[1] <- S
  ObsS <- array(dim=nyears+1)
  randS <- rnorm(1)
  ObsS[1] <- max(spawners[1]*exp(SCV*randS),0.) 
  
  recruits <- array(dim=nyears)
  catch <- array(dim=nyears)
  zold <- 0
  
  
  #### actual simulation with goal adjustments #############################
  sims <- paste(outname,"Sims.txt")
  sims <- gsub(" ", "", sims, fixed = TRUE)
  sumsims <- paste(outname,"Sum.txt")
  sumsims <- gsub(" ", "", sumsims, fixed = TRUE)
  write(c("iter   i     a     b      Smsy      Emin"),file=sims)
  write(c("iter   avgS    lowS%    avgR     lowR      avgC      noC"),file=sumsims)
  for (iter in 1:niters) {
    #  browser()
    cat(iter,"\n") 
    if(SRmodel=="B") {
      trueSmsy <- BSmsy(alpha,beta)
    } else {
      trueSmsy <- RSmsyExact(c(alpha,beta))$Smsy
    }
    spawners[1] <- trueSmsy*runif(1,min=0.5,max=1.5)
    S <- spawners[1]
    ObsS <- array(dim=nyears+1)
    randS <- rnorm(1)
    ObsS[1] <- max(spawners[1]*exp(SCV*randS),0.) 
    recruits <- array(dim=nyears)
    catch <- array(dim=nyears)
    zold <- 0
    if(management=="precise") {
      Emin <- trueSmsy #constant escapement goal policy
      hrate <- 1
    } else {
      Emin <- 0.8*trueSmsy #constant escapement goal policy
      hrate <- 0.8
    }
    if(harvHist=="Overh") {
      Emin=Emin/2
    } else if (harvHist=="Underh") {
      Emin=Emin*1.5
    }
    
    for (i in 1:nyears) {
      #if((i>=leapyear)&(i<(leapyear+leapduration))) {zold <- leapamt}
      alpha_yr <- ifelse((i>=leapyear)&(i<(leapyear+leapduration)),alpha+leapa,alpha) 
      beta_yr <- ifelse((i>=leapyear)&(i<(leapyear+leapduration)),beta*leapb,beta) 
      if(SRmodel=="B") {
        simout <- simBH(alpha_yr,beta_yr,phi,RCV,SCV,Emin,hrate,CVE,spawners[i])
      } else {
        simout <- simRicker(alpha_yr,beta_yr,phi,RCV,SCV,Emin,hrate,CVE,spawners[i])
      }
      #    cat(i,simout,"\n")
      recruits[i] <- simout[1]
      spawners[i+1] <- simout[2]
      S <- simout[2]
      ObsS[i+1] <- simout[3]
      catch[i] <- simout[4]
      zold <- simout[5]
      #    write(c(iter,i,spawners[i],recruits[i]),file="junk.txt",append=TRUE)
      
      if(NewGoal[i]) {
        # fit simulated data
        vars <- list(n=i-1,s=ObsS[1:(i-1)],r=recruits[1:(i-1)])
        pars <- c(alpha,beta)
        if(SRmodel=="B") {
          z <- optim(pars,Bfnorm,data=vars)
        } else {
          z <- optim(pars,Rfnorm,data=vars)
        }
        afit <- abs(z$par[1])
        bfit <- abs(z$par[2])
        if(SRmodel=="B") {
          sfit <- BSmsy(afit,bfit)
        } else {
          sfit <- RSmsyExact(c(afit,bfit))$Smsy
        }
        #update escapement goal
#        rediff <- abs(sfit-Emin)/Emin
#        diff <- sign(sfit-Emin)
#        if(rediff>max_goal_change) { # prevent large changes in goal
#          Emin <- Emin*(1+diff*max_goal_change)
#        } else { # acceptable amount of change
#          Emin <- sfit
#        }
#        if(management!="precise") Emin <- 0.8*Emin
        
        z$par
        #cat(iter,i,abs(z$par),sfit,Emin,"\n")
        write(c(iter,i,abs(z$par),sfit,Emin),ncolumns=6,file=sims,append=TRUE)
        
        #simdat <- data.frame(spawners=vars$s,recruits=vars$r,years=seq(from=1,to=i-1))
        #pp <- ggplot(data=simdat)+geom_line(aes(years,spawners))
        #pp <- pp+geom_line(aes(years,recruits),color="red")+ylim(0,NA)
        #pp
        
        
      } #new goal - refit data
      
    } # for i loop
    avgS <- mean(spawners)
    lowS <- min(spawners)/trueSmsy
    avgR <- mean(spawners)
    lowR <- min(spawners)
    avgC <- mean(catch)
    noC <- sum(catch==0)
    write(c(iter,avgS,lowS,avgR,lowR,avgC,noC),ncolumns=7,file=sumsims,append=TRUE)
    
    
  } #### iter - main simulation
  
} #dummy loop


result <- as.data.frame(read.table(sims,header=TRUE))

mm <- aggregate(result,by=list(result$i),median)
#plot(mm[,1],mm[,6])
df <- data.frame(year=mm[,1],avg_goal=mm[,5])
aggregate(result,by=list(result$i),mean)
xx <- aggregate(result,by=list(result$i),sd)
ss <- xx/sqrt(niters)
ss

df <- data.frame(year=as.ordered(result[,2]),egoal=result[,5],alpha=result[,3],beta=result[,4])
pp <- ggplot(data=df,aes(y=egoal,x=year))+geom_boxplot()+ylim(0,3*median(df$egoal))
pp <- pp+geom_abline(slope=0,intercept=trueSmsy,color="red",size=1.5)
if(leapb != 1) pp <- pp+geom_abline(slope=0,intercept=leapSmsy,color="blue",size=1.5)
pp <- pp+theme(text = element_text(size = 20)) 
pp
fname <- paste(c(outname,"_escgoal.jpg"),collapse="")
fname
dev.copy(jpeg,fname)
dev.off()
pp <- ggplot(data=df,aes(y=beta,x=year))+geom_boxplot()+ylim(0,3*median(df$beta))
pp <- pp+geom_abline(slope=0,intercept=beta,color="red",size=1.5)
if(leapb != 1) pp <- pp+geom_abline(slope=0,intercept=beta*leapb,color="blue",size=1.5)
pp <- pp+theme(text = element_text(size = 20)) 
pp
fname <- paste(c(outname,"_beta.jpg"),collapse="")
fname
dev.copy(jpeg,fname)
dev.off()
pp <- ggplot(data=df,aes(y=alpha,x=year))+geom_boxplot()+ylim(0,3*median(df$alpha))
pp <- pp+geom_abline(slope=0,intercept=alpha,color="red",size=1.5)
if(leapb != 1) pp <- pp+geom_abline(slope=0,intercept=alpha+leapa,color="blue",size=1.5)
pp <- pp+theme(text = element_text(size = 20)) 
pp
fname <- paste(c(outname,"_alpha.jpg"),collapse="")
fname
dev.copy(jpeg,fname)
dev.off()

