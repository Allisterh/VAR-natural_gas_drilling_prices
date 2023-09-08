##############################################################################
### Script: Bayes 4x4 VAR
### Author: Valentin Winkler
### Date:   08/2022
### Description:  Read in Data for 4x4 gas-market VAR with exogenos variables,
###               manipulate data, get posterior distribution of reduced form
###               coefficients and after that draw M representations of reduced
###               form parameters and for each representation draw K orthogonal
###               matrices Q satisfying the zero restrictions and check if sign 
###               restrictions are met.
###               Afterwards construct joint posterior mode and 68% credibility
###               band like in Kilian and Inoue (2013) for structural IRFs.
##############################################################################
### Inputs:   Start and End Date of estimation window, draws of reduced form 
###           parameters and orthogonal matrices
### optional:
#   a34       restriction of ipd reaction to gas demand shock on impact (default 0.004)
#   ela       restriction of price elasticity of gas supply on impact (default 0.065)
#   tren      logical, indicator if broken trend is included (default TRUE)
###
### Outputs:  IRFs, model parameters and posterior density value of admissible
###           draws
##############################################################################

############## Packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(vars)
library(ggpubr)
library(reshape2)
library(scales)
source("bayes_rotations_4_4.R")

bayes_4_4_var <- function(star, en, M, K, a34 = 0.004, ela = 0.065, tren = TRUE){
  ##############################################################################
  ############## Parameters ####################################################
  ##############################################################################
  
  # specification and simulation parameters
  p <- 6 # lag of VAR model
  h <- 12 # horizon of IRFs to be computed
  
  # number of dummy variables 
  d <- 3
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = M, initial = 0, style = 3)
  
  ### prior parameters are specified below
  
  ##############################################################################
  ############## Read in and manipulate data ###################################
  ##############################################################################
  
  D <- read.csv("Data_VAR.csv", sep = ";", dec = ",")
  
  #### Manipulate Data ####
  D$Date <- as.Date(D$Date, format = "%d.%m.%Y") # date format
  
  # make longer Henry Hub Time Series
  D$gasprice_norm[411:510] <- D$gasprice_norm2[411:510]
  
  # compute real gas price
  D <- D %>%
    mutate(gasprice = gasprice_norm / cpi)
  
  # for heating/cooling degree days: get average of last 5 years
  D$heat_rollm <- NA
  D$cool_rollm <- NA
  
  for (i in (5*12 +1):dim(D)[1]) {
    D$heat_rollm[i] <- mean(D$heating_days[i-5*12],
                            D$heating_days[i-4*12],
                            D$heating_days[i-3*12],
                            D$heating_days[i-2*12],
                            D$heating_days[i-1*12])
    
    D$cool_rollm[i] <- mean(D$cooling_days[i-5*12],
                            D$cooling_days[i-4*12],
                            D$cooling_days[i-3*12],
                            D$cooling_days[i-2*12],
                            D$cooling_days[i-1*12])
  }
  rm(i)
  
  # get deviation of heating/cooling days from 5Y average
  D <- D %>%
    mutate(heat_dev = heating_days - heat_rollm,
           cool_dev = cooling_days - cool_rollm)
  
  # transform rigs, gas production, ind. prod. as well as real
  # price of gas to log levels
  D <- D %>%
    mutate(rig = log(rig_count),
           gpd = log(gasprod),
           ipd = log(ip),
           rpg = log(gasprice))
  
  # data frame of endogenous and exogenous variables + dataset of dummies
  end <- subset(D, Date > star & Date < en) %>% dplyr::select(rig, gpd, ipd, rpg)
  exo <- subset(D, Date > star & Date < en) %>% dplyr::select(heat_dev, cool_dev)
  dum <- subset(D, Date > star & Date < en) %>% dplyr::select(Date)
  
  # create dummies
  dum$dum1 <- as.numeric(dum$Date == as.Date("2005-09-15")) # Katrina dummy
  dum$dum2 <- as.numeric(dum$Date == as.Date("2008-09-15")) # Gustav dummy
  dum$dum3 <- as.numeric(dum$Date >= as.Date("2005-05-15"))
  dum$dum3 <- cumsum(dum$dum3)
  
  if(tren == FALSE){ # set trend dummy to zero (will not be included in regressor matrix automatically)
    dum$dum3 <- dum$dum3*0
  }
  
  # add p lags of endogenous variables to data
  name_tmp <- names(end)
  
  for (i in 1:p) {
    end <- data.frame(end, lag(end$rig, n = i), lag(end$gpd, n = i), lag(end$ipd, n = i), lag(end$rpg, n = i))
    
    # new variable names
    txt_rig <- paste("rig_l",i,sep = "")
    txt_gpd <- paste("gpd_l",i,sep = "")
    txt_ipd <- paste("ipd_l",i,sep = "")
    txt_rpg <- paste("rpg_l",i,sep = "")
    
    name_tmp <- c(name_tmp, txt_rig, txt_gpd, txt_ipd, txt_rpg)
  }
  
  # rename data frame
  names(end) <- name_tmp
  rm(name_tmp, txt_rig, txt_gpd, txt_ipd, txt_rpg, i)
  
  # remove first p observations
  end <- end[-c(1:p),]
  exo <- exo[-c(1:p),]
  dum <- dum[-c(1:p),]
  
  ### create list to store admissible IRFs and corresponding parameters in
  resu <- list(
    irig = list(
      rrig = data.frame(Horizon = 0:h),
      rgpd = data.frame(Horizon = 0:h),
      ripd = data.frame(Horizon = 0:h),
      rrpg = data.frame(Horizon = 0:h)
    ),
    igpd = list(
      rrig = data.frame(Horizon = 0:h),
      rgpd = data.frame(Horizon = 0:h),
      ripd = data.frame(Horizon = 0:h),
      rrpg = data.frame(Horizon = 0:h)
    ),
    iipd = list(
      rrig = data.frame(Horizon = 0:h),
      rgpd = data.frame(Horizon = 0:h),
      ripd = data.frame(Horizon = 0:h),
      rrpg = data.frame(Horizon = 0:h)
    ),
    irpg = list(
      rrig = data.frame(Horizon = 0:h),
      rgpd = data.frame(Horizon = 0:h),
      ripd = data.frame(Horizon = 0:h),
      rrpg = data.frame(Horizon = 0:h)
    ),
    BET = list(), # reduced form slope parameters
    C = list(), # Chol(SIG)
    Q = list(), # orthogonal matrices
    poster_f = c() # posterior density values
  )
  
  # Get length of time horizon
  T <- dim(end)[1]
  T_y <- ceiling(T/12) # get number of years ronded to next integer
  
  # Get 4 x T matrix of dependent variables
  Y <- t(as.matrix(end[,1:4]))
  
  # Construct (1 (Intercept) + 4p (Endog.) + 2 (Exog.) + 11 (Season) + d (dummies)) x T
  # regressor matrix 
  Z <- data.frame(Inter = rep(1,T), # intercept
                  end[, -c(1:4)], # lagged endogenous variables
                  exo) 
  
  # remove non-transformed data
  rm(D, end, exo)
  
  # Create 11 seasonal dummies
  Z$S1 <- rep(c(1,0,0,0,0,0,0,0,0,0,0,0),T_y)[1:T]
  Z$S2 <- rep(c(0,1,0,0,0,0,0,0,0,0,0,0),T_y)[1:T]
  Z$S3 <- rep(c(0,0,1,0,0,0,0,0,0,0,0,0),T_y)[1:T]
  Z$S4 <- rep(c(0,0,0,1,0,0,0,0,0,0,0,0),T_y)[1:T]
  Z$S5 <- rep(c(0,0,0,0,1,0,0,0,0,0,0,0),T_y)[1:T]
  Z$S6 <- rep(c(0,0,0,0,0,1,0,0,0,0,0,0),T_y)[1:T]
  Z$S7 <- rep(c(0,0,0,0,0,0,1,0,0,0,0,0),T_y)[1:T]
  Z$S8 <- rep(c(0,0,0,0,0,0,0,1,0,0,0,0),T_y)[1:T]
  Z$S9 <- rep(c(0,0,0,0,0,0,0,0,1,0,0,0),T_y)[1:T]
  Z$S10 <- rep(c(0,0,0,0,0,0,0,0,0,1,0,0),T_y)[1:T]
  Z$S11 <- rep(c(0,0,0,0,0,0,0,0,0,0,1,0),T_y)[1:T]
  
  d2 <- 0 # number of dummy variables really included
  
  for (j in 1:d) {
    if(var(dum[ ,1+j]) != 0){ # if dummy column not column of zeros/ones...
      Z <- data.frame(Z, dum[ ,1+j]) # ...include in regressor matrix
      d2 <- d2 + 1
    }
  }
  
  Z <- t(as.matrix(Z))
  
  ##############################################################################
  ############## Compute posterior distribution ################################
  ##############################################################################
  
  # prior parameters to get LS point estimator
  Betbar0 <- matrix(0, nrow = 4, ncol = (1 + 4*p + 2 + 11 + d2))
  nu0 <- 0
  V1 <- matrix(0, nrow = (1 + 4*p + 2 + 11 + d2), ncol = (1 + 4*p + 2 + 11 + d2))
  S0 <- matrix(0, nrow = 4, ncol = 4)
  
  # Compute LS estimators of regression coefficients
  BET <- Y %*% t(Z) %*% solve(Z %*% t(Z))
  
  # coefficients corresponding to autoregressive terms
  A <- BET[,2:(p*4 + 1)]
  
  # compute residual matrix
  U <- Y - BET %*% Z
  
  # compute estimator for innovation covariance matrix
  SIG <- ( U %*% t(U) )/T
  
  ####### Posterior parameters
  nuT <- T + nu0
  VT <- V1 + Z %*% t(Z)
  BetbarT <- (Betbar0 %*% V1 + Y %*% t(Z)) %*% solve(V1 + Z %*% t(Z))
  ST <- T * SIG + S0 + BET %*% Z %*% t(Z) %*% t(BET) + Betbar0 %*% V1 %*% t(Betbar0) - BetbarT %*% (V1 + Z %*% t(Z)) %*% t(BetbarT)
  betbarT <- matrix(BetbarT, ncol = 1)
  
  ##############################################################################
  ############## Draw from posterior distribution ##############################
  ##############################################################################
  for (i in 1:M) {
    # update progress bar
    setTxtProgressBar(pb,i)
    
    # compute sigma
    SIG_draw <- matrix(0,nrow=4, ncol=4)
    for (j in 1:nuT) {
      xi <- t(chol(solve(ST))) %*% matrix(rnorm(4), nrow =4)
      SIG_draw <- SIG_draw + xi %*% t(xi)
    }
    SIG_draw <- solve(SIG_draw)
    
    # VC-Matrix of vec(BET)
    SIG_bet <- solve(VT) %x% SIG_draw
    
    bet_hat <- betbarT + t(chol(SIG_bet)) %*% matrix(rnorm(4*(1 + 4*p + 2 + 11 + d2)), ncol = 1)
    BET_hat <- matrix(bet_hat, nrow = 4)
    
    # Generate K realizations of orthogonal matrices Q and check if they fulfill
    # sign restrictions
    resu_tmp <- bayes_rotations_4_4(K, h, p, BET_hat, SIG_draw, betbarT, VT, nuT, ST, a34, ela)
    
    # Append results for current posterior draw to final result matrix
    resu$irig$rrig <- data.frame(resu$irig$rrig, resu_tmp$irig$rrig[ ,-1])
    resu$irig$rgpd <- data.frame(resu$irig$rgpd, resu_tmp$irig$rgpd[ ,-1])
    resu$irig$ripd <- data.frame(resu$irig$ripd, resu_tmp$irig$ripd[ ,-1])
    resu$irig$rrpg <- data.frame(resu$irig$rrpg, resu_tmp$irig$rrpg[ ,-1])
    resu$igpd$rrig <- data.frame(resu$igpd$rrig, resu_tmp$igpd$rrig[ ,-1])
    resu$igpd$rgpd <- data.frame(resu$igpd$rgpd, resu_tmp$igpd$rgpd[ ,-1])
    resu$igpd$ripd <- data.frame(resu$igpd$ripd, resu_tmp$igpd$ripd[ ,-1])
    resu$igpd$rrpg <- data.frame(resu$igpd$rrpg, resu_tmp$igpd$rrpg[ ,-1])
    resu$iipd$rrig <- data.frame(resu$iipd$rrig, resu_tmp$iipd$rrig[ ,-1])
    resu$iipd$rgpd <- data.frame(resu$iipd$rgpd, resu_tmp$iipd$rgpd[ ,-1])
    resu$iipd$ripd <- data.frame(resu$iipd$ripd, resu_tmp$iipd$ripd[ ,-1])
    resu$iipd$rrpg <- data.frame(resu$iipd$rrpg, resu_tmp$iipd$rrpg[ ,-1])
    resu$irpg$rrig <- data.frame(resu$irpg$rrig, resu_tmp$irpg$rrig[ ,-1])
    resu$irpg$rgpd <- data.frame(resu$irpg$rgpd, resu_tmp$irpg$rgpd[ ,-1])
    resu$irpg$ripd <- data.frame(resu$irpg$ripd, resu_tmp$irpg$ripd[ ,-1])
    resu$irpg$rrpg <- data.frame(resu$irpg$rrpg, resu_tmp$irpg$rrpg[ ,-1])
    
    resu$BET <- append(resu$BET, resu_tmp$BET)
    resu$C <- append(resu$C, resu_tmp$C)
    resu$Q <- append(resu$Q, resu_tmp$Q)
    resu$poster_f <- c(resu$poster_f, resu_tmp$poster_f)
  }
  
  # Append other model parameters to result list
  resu <- append(resu, list(LS_model = list(Z = Z, Y = Y, SIG = SIG, BET = BET)))
  
  close(pb) # close progress bar
  return(resu)
}
