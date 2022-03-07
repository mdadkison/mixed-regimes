###################################################################
simRicker <- function(alpha,beta,phi,RCV,SCV,Emin,hrate,ECV,zold) { 
  # NO AGE STRUCTURE - ASSUMES ONE AGE AT RETURN
  #  projects population forward one generation using a Ricker model
  #  with log-normal, AR(1) process error in recruits and
  #  log-normal observation error in spawners and
  #  log-normal error in achieving escapement target
  #  for equation R <- (maxRpS*S)/(1.+maxRpS*S/capacity)
  #  Note - all log-normal error terms adjusted for upwards bias
  
  #inputs:
  #  alpha = ln(R/S) at low stock size #######NOTE alpha != maxRpS ########
  #  beta = unfished equilibrium stock size = beta
  #  phi = AR(1) autocorrelation
  #  RCV = std dev of log-normal error in recruitment
  #  SCV = std dev of log-normal error in observing spawners
  #  Emin = target minimum escapement
  #  hrate = target harvest rate on recruits above minimum escapement
  #  ECV = std dev of log-normal implementation error in achieving desired spawners
  #  S = number of spawners producing this year's recruitment 
  #  zold = previous recruitment deviation
  #
  #outputs function value as a vector:
  #  (recruits,spawners,ObsS,catch,zold)
  #  recruits = true number of offspring produced by S
  #  spawners = true number of spawners after recruits reduced by catch
  #  ObsS = observed value for number of spawners
  #  catch = true value for catch from the recruits
  #  zold = deviation of recruitment from the Ricker expectation
  ########################################################
  #browser()
  # simulate data
  #recruits
  #alpha <- maxRpS/exp(1) #NOTE alpha != maxRpS
  expect <- S*exp(alpha*(1-S/beta))
  randR <- rnorm(1)
  z <- ((1-phi)*randR+phi*zold)/((1-phi)*(1-phi)+phi*phi)
  #cat(i,randR,z,zold,"\n")
  recruits <- max(expect*exp(RCV*z-0.5*RCV*RCV),0.)
  #cat(i,randR,z,expect,recruits,"\n")
  zold <- log(recruits)-log(expect) #residual for next year's autocorrelation
  #spawners
  espawn <- max(Emin + (1-hrate)*(recruits-Emin),0) #expected spawners from HCR
  randE <- rnorm(1)
  spawners <- min(espawn*exp(ECV*randE-0.5*ECV*ECV),recruits)
  #observed spawners
  randS <- rnorm(1)
  ObsS <- max(spawners*exp(SCV*randS-0.5*SCV*SCV),0.) #observed escapement
  #catch
  catch <- max(recruits-spawners,0.)
  #  cat(i,randR,z,spawners,recruits,catch,"\n")
  #  }
  out <- c(recruits,spawners,ObsS,catch,zold)
  return(out)
} ########## function simRicker #######################

###################################################################
simBH <- function(alpha,beta,phi,RCV,SCV,Emin,hrate,ECV,S,zold) { 
  # NO AGE STRUCTURE - ASSUMES ONE AGE AT RETURN
  #  projects population forward one generation using a Ricker model
  #  with log-normal, AR(1) process error in recruits and
  #  log-normal observation error in spawners and
  #  log-normal error in achieving escapement target
  #  for equation R <- (maxRpS*S)/(1.+maxRpS*S/capacity)
  #  Note - all log-normal error terms adjusted for upwards bias
  
  #inputs:
  #  capacity = unfished equilibrium stock size ########NOTE - beta != capacity ################
  #  maxRpS = R/S at low stock size
  #  phi = AR(1) autocorrelation
  #  RCV = std dev of log-normal error in recruitment
  #  SCV = std dev of log-normal error in observing spawners
  #  Emin = target minimum escapement
  #  hrate = target harvest rate on recruits above minimum escapement
  #  ECV = std dev of log-normal implementation error in achieving desired spawners
  #  S = number of spawners producing this year's recruitment 
  #  zold = previous recruitment deviation
  #
  #outputs function value as a vector:
  #  (recruits,spawners,ObsS,catch,zold)
  #  recruits = true number of offspring produced by S
  #  spawners = true number of spawners after recruits reduced by catch
  #  ObsS = observed value for number of spawners
  #  catch = true value for catch from the recruits
  #  zold = deviation of recruitment from the Ricker expectation
  ########################################################
  #browser()
  # simulate data
  #recruits
  #beta <- capacity*maxRpS/(maxRpS-1) #NOTE - beta != capacity
  expect <- (alpha*S)/(1.+alpha*S/beta)
  randR <- rnorm(1)
  z <- ((1-phi)*randR+phi*zold)/((1-phi)*(1-phi)+phi*phi)
  #cat(i,randR,z,zold,"\n")
  recruits <- max(expect*exp(RCV*z-0.5*RCV*RCV),0.)
  #cat(i,randR,z,expect,recruits,"\n")
  zold <- log(recruits)-log(expect) #residual for next year's autocorrelation
  #spawners
  espawn <- max(Emin + (1-hrate)*(recruits-Emin),0) #expected spawners from HCR
  randE <- rnorm(1)
  spawners <- min(espawn*exp(ECV*randE-0.5*ECV*ECV),recruits)
  #observed spawners
  randS <- rnorm(1)
  ObsS <- max(spawners*exp(SCV*randS-0.5*SCV*SCV),0.) #observed escapement
  #catch
  catch <- max(recruits-spawners,0.)
  #  cat(i,randR,z,spawners,recruits,catch,"\n")
  #  }
  out <- c(recruits,spawners,ObsS,catch,zold)
  return(out)
} ########## function simRicker #######################


# @@@@@@@ function returns SSQ from log-transformed data @@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@ given proposed parameter values @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Bflognorm <- function(pars,data) {
  #calculates ssq for Ricker fit using ln(R/S)
  #browser()
  alpha <- abs(pars[1])
  beda <- abs(pars[2])
  n <- data$n   
  S <- data$s
  R <- data$r
  ssq <- 0.
  for (i in 1:n) {
    predict <- log((alpha*S[i])/(1.+alpha*S[i]/beda))
    ssq <- ssq + (log(R[i])-predict)^2
  } #i
  return(ssq)
} 
Rflognorm <- function(pars,data) {
  #calculates ssq for Ricker fit using ln(R/S)
  #browser()
  alpha <- abs(pars[1])
  beta <- abs(pars[2])
  n <- data$n   
  S <- data$s
  R <- data$r
  ssq <- 0.
  for (i in 1:n) {
    predict <- log(S[i])+alpha*(1-S[i]/beta)
    ssq <- ssq + (log(R[i])-predict)^2
  } #i
  return(ssq)
} 

#@@@@@@@@@@@@@@@@@ function fr   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

# @@@@@@@ function returns SSQ from normally-distributed data @@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@ given proposed parameter values @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Bfnorm <- function(pars,data) {
  #calculates ssq for Ricker fit using ln(R/S)
  #browser()
  alpha <- abs(pars[1])
  beda <- abs(pars[2])
  n <- data$n   
  S <- data$s
  R <- data$r
  ssq <- 0.
  for (i in 1:n) {
    predict <- (alpha*S[i])/(1.+alpha*S[i]/beda)
    ssq <- ssq + (R[i]-predict)^2
  } #i
  return(ssq)
} 
Rfnorm <- function(pars,data) {
  #calculates ssq for Ricker fit using ln(R/S)
  #browser()
  alpha <- abs(pars[1])
  beta <- abs(pars[2])
  n <- data$n   
  S <- data$s
  R <- data$r
  ssq <- 0.
  for (i in 1:n) {
    predict <- (S[i])*exp(alpha*(1-S[i]/beta))
    ssq <- ssq + (R[i]-predict)^2
  } #i
  return(ssq)
} 

#@@@@@@@@@@@@@@@@@ function fr   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

BSmsy <- function(alpha,beta) {
  #NOTE - beta <> capacity
  SS <- beta*(sqrt(1./alpha)-1./alpha)
  return(SS)
}

RSmsy <- function(alpha,beta) {
  # R = S*exp(alpha*(1-S/beta))
  # Hilborn approximation
  SS <- beta*(0.5-0.07*alpha)
  return(SS)
}

RSmsyExact <- function(pars) {
  # R = S*exp(alpha*(1-S/beta))
  #uses optimise to get the exact Smsy rather than the Hilborn approximation
  Ryield <- function(S,pars) {
    # gives the surplus Ricker yield for a given spawner abundance
    a <- pars[1]
    b <- pars[2]
    Y <- S*exp(a*(1-S/b))-S
    return(Y)
  } #internal function Ryield
  zz <- optimise(Ryield,lower=0,upper=pars[2],pars=pars,maximum=TRUE)
  outlist <- list(MSY=zz$objective,Smsy=zz$maximum
                  ,hmsy=zz$objective/(zz$objective+zz$maximum))
  return(outlist)
} #SmsyExact


Bhmsy <- function(alpha,beta) {
  #NOTE - beta <> capacity
  HH <- 1.-sqrt(1./alpha)
  return(HH)
}  

Rhmsy <- function(alpha,beta) {
  # R = S*exp(alpha*(1-S/beta))
  # Hilborn approximation
  HH <- 0.5*alpha-0.07*alpha*alpha
  return(HH)
}  


