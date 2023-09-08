##############################################################################
### Script: Bayes rotations 3x3
### Author: Valentin Winkler
### Date:   08/2022
### Description:  Get estimates of reduced form VAR-coefficients for gas market
###               model with 4 endogenous variables
###               Draw K orthogonal matrices Q satisfying zero restrictions and 
###               compute IRFs of horizon h if B_0^-1 = chol(SIGMA)Q and IRFs 
###               satisfy sign restrictions
###               Also, compute the posterior density of IRFs
##############################################################################

######## Input
# K       Number of draws of orthogonal matrix Q
# h       Horizon of IRFs to be computed (must be >= 6)
# p       Lags of the VAR model
# BET     4 x (1 + 4p + 2 + 11) Matrix of estimated, reduced form slope parameters
# SIG     4 x 4 Matrix, estimated covariance-matrix of reduced form errors
# betbarT posterior mean of vec(BET)
# VT      V1 + Z*Z'
# nuT     df for inverse-Wishart distribution
# ST      Scale matrix for inverse-Wishart distribution
#
### optional:
# a34     restriction of ipd reaction to gas demand shock on impact (default 0.004)
# ela     restriction of price elasticity of gas supply on impact (default 0.065)

######## Output
# List 'resu' with...
#     Lists irig/igpd/iipd/irpg that contain responses of all variables to structural
#                           shocks to rig/gpd/ipd/rpg
#     List A that contains matrices of estimated slope parameters
#     List C that contains chol(SIG) for every draw
#     List Q that contains every admissible Q
#     Vector poster_f with posterior values of draw

##############################################################################
### Exclusion restrictions and sign/elasticity restrictions imposed
### Variable Ordering:  rig (active gas drilling rigs), gpd (gas production), 
###                     ipd (industrial production), rpg (real price of gas)

# Impact:
#     * 0 0 0
#     * - + +
#     * - + -
#     * + + +

# Elasticity constraints
# -  The Impact price elasticity of gas production is bound to 0.065
#    This is motivated by Mason and Roberts (2018) who estimate that the impact
#    price elasticity of gas production of existing wells is below 0.03 on average
#    and slightly below 0.065 for the 25% most productive well.
#    Since drilling of new wells takes approximately three months to even get 
#    started (Kellogg 2014) and after that at least another month to start production,
#    short run gas supply is bounded by the elasticity of existing wells.
# -  An one standard deviation gas demand shock is allowed to decrease industrial
#    production at most by 0.004 (0.008 is the standard deviation of the reduced
#    form error of ipd)

##############################################################################


##### Load Packages #####
library(dplyr)
library(tidyr)
source("irfpdf_blockrec.R")

bayes_rotations_4_4 <- function(K, h, p, BET, SIG, betbarT, VT, nuT, ST, a34 = 0.004, ela = 0.065){
  
  # cholesky decomposition of SIGMA
  C <- t(chol(SIG))
  
  # get AR-coefficients
  A <- BET[,2:(p*4 + 1)]
  
  # Submatrix for rotations
  Csub <- C[-1,-1]
  
  # for convenience: store rows of C as row vectors
  c1 <- Csub[1,]
  c2 <- Csub[2,]
  c3 <- Csub[3,]
  
  # construct companion matrix for VAR model
  AA <- matrix(0, 4*p, 4*p)
  AA[1:4,] <- A
  AA[-c(1:4),1:((p-1)*4)] <- diag((p-1)*4)
  
  # J-Matrix to get IRF as J (AA)^i J'
  J <- matrix(0,4,4*p)
  J[1:4,1:4] <- diag(4)
  
  ### Compute reduced form IRFs
  RIRF <- list(C, J %*% AA %*% t(J) %*% C)
  tmp <- AA %*% AA
  
  for (i in 2:h) {
    RIRF[[1+i]] <- (J %*% tmp %*% t(J) %*% C)
    tmp <- tmp %*% AA
  }
  rm(tmp,i)
  
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
    poster_f = c() # posterior value of draw
  )
  
  # go through K random draws for B0^-1
  for (i in 1:K) {
    
    # random 3x3-Matrix with IID N(0,1) entries
    W <- matrix(rnorm(9),3)
    QR <- qr(W)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    
    # make sure that diagonal entries of R are positive, otherwise flip signs
    for (j in 1:3) {
      if(R[j,j] < 0)
        Q[,j] <- -Q[,j]
    }
    rm(j, R, QR, W)
    
    # for convenience: store columns of Q as column vectors
    q1 <- Q[,1]
    q2 <- Q[,2]
    q3 <- Q[,3]
    
    #####################################################################
    ### Check if impact sign restrictions are met
    ### if not: throw away draw
    ### for efficiency: check if columns can be switched to fulfill sign restrictions
    #####################################################################
    
    ###### Case 1: Column 1 ok
    if((c1 %*% q1) < 0 & (c2 %*% q1) < 0 & (c3 %*% q1) > 0){ # 1st column ok
      
      if((c1 %*% q2) > 0 & (c2 %*% q2) > 0 & (c3 %*% q2) > 0){ # 2nd column ok
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q2)) > 0 & (c2 %*% (-q2)) > 0 & (c3 %*% (-q2)) > 0){ # -q2 ok
        
        q2 <- -q2 # switch sign of column
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% q3) > 0 & (c2 %*% q3) > 0 & (c3 %*% q3) > 0){ # q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) > 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- -q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else{
        next # throw away draw
      }
    }
    
    ###### Case 2: -q1 ok
    else if((c1 %*% (-q1)) < 0 & (c2 %*% (-q1)) < 0 & (c3 %*% (-q1)) > 0){ # -q1 ok
      
      q1 <- -q1 # change sign of column
      
      if((c1 %*% q2) > 0 & (c2 %*% q2) > 0 & (c3 %*% q2) > 0){ # 2nd column ok
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q2)) > 0 & (c2 %*% (-q2)) > 0 & (c3 %*% (-q2)) > 0){ # -q2 ok
        
        q2 <- -q2 # switch sign of column
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% q3) > 0 & (c2 %*% q3) > 0 & (c3 %*% q3) > 0){ # q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) > 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- -q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else{
        next # throw away draw
      }
    }
    
    ###### Case 3: q2 ok as q1
    else if((c1 %*% (q2)) < 0 & (c2 %*% (q2)) < 0 & (c3 %*% (q2)) > 0){ # q2 ok as q1
      
      tmp <- q1 # change columns
      q1 <- q2
      q2 <- tmp
      
      if((c1 %*% q2) > 0 & (c2 %*% q2) > 0 & (c3 %*% q2) > 0){ # 2nd column ok
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q2)) > 0 & (c2 %*% (-q2)) > 0 & (c3 %*% (-q2)) > 0){ # -q2 ok
        
        q2 <- -q2 # switch sign of column
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% q3) > 0 & (c2 %*% q3) > 0 & (c3 %*% q3) > 0){ # q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) > 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- -q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else{
        next # throw away draw
      }
    }
    
    ###### Case 4: -q2 ok as q1
    else if((c1 %*% (-q2)) < 0 & (c2 %*% (-q2)) < 0 & (c3 %*% (-q2)) > 0){ # -q2 ok as q1
      
      tmp <- q1 # change columns
      q1 <- -q2
      q2 <- tmp
      
      if((c1 %*% q2) > 0 & (c2 %*% q2) > 0 & (c3 %*% q2) > 0){ # 2nd column ok
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q2)) > 0 & (c2 %*% (-q2)) > 0 & (c3 %*% (-q2)) > 0){ # -q2 ok
        
        q2 <- -q2 # switch sign of column
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% q3) > 0 & (c2 %*% q3) > 0 & (c3 %*% q3) > 0){ # q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) > 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- -q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else{
        next # throw away draw
      }
    }
    
    ###### Case 5: q3 ok as q1
    else if((c1 %*% (q3)) < 0 & (c2 %*% (q3)) < 0 & (c3 %*% (q3)) > 0){ # q3 ok as q1
      
      tmp <- q1 # change columns
      q1 <- q3
      q3 <- tmp
      
      if((c1 %*% q2) > 0 & (c2 %*% q2) > 0 & (c3 %*% q2) > 0){ # 2nd column ok
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q2)) > 0 & (c2 %*% (-q2)) > 0 & (c3 %*% (-q2)) > 0){ # -q2 ok
        
        q2 <- -q2 # switch sign of column
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% q3) > 0 & (c2 %*% q3) > 0 & (c3 %*% q3) > 0){ # q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) > 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- -q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else{
        next # throw away draw
      }
    }
    
    ###### Case 6: -q3 ok as q1
    else if((c1 %*% (-q3)) < 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q1
      
      tmp <- q1 # change columns
      q1 <- -q3
      q3 <- tmp
      
      if((c1 %*% q2) > 0 & (c2 %*% q2) > 0 & (c3 %*% q2) > 0){ # 2nd column ok
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q2)) > 0 & (c2 %*% (-q2)) > 0 & (c3 %*% (-q2)) > 0){ # -q2 ok
        
        q2 <- -q2 # switch sign of column
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% q3) > 0 & (c2 %*% q3) > 0 & (c3 %*% q3) > 0){ # q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) > 0 & (c3 %*% (-q3)) > 0){ # -q3 ok as q2
        tmp <- q2 # Switch columns
        q2 <- -q3
        q3 <- tmp
        
        if((c1 %*% q3) > 0 & (c2 %*% q3) < 0 & (c3 %*% q3) > 0){ # 3rd column ok
          q3 <- q3 # keep draw
        }
        else if((c1 %*% (-q3)) > 0 & (c2 %*% (-q3)) < 0 & (c3 %*% (-q3)) > 0){ # -q3 ok
          q3 <- -q3 # keep draw
        }
        else{
          next # throw away draw
        }
      }
      
      else{
        next # throw away draw
      }
    }
    
    ##### Case 7: No column ok as q1
    else{
      next # throw away draw
    }
    
    #########################################################################
    ######################### Check dynamic sign restrictions ###############
    #########################################################################
    
    # Construct full, 4x4 orthogonal matrix
    Q <- diag(4)
    Qsub <- cbind(q1, q2, q3)
    Q[-1,-1] <- Qsub
    
    ######### DYNAMIC RESTRICTION 1
    ######### In first 6 months:  economic activity shock should raise ind. prod.
    
    # economic activity
    ind <- 0
    
    #for (j in 0:6) {
    #  if((RIRF[[j+1]][3, ] %*% Q[,3] < 0)){
    #    ind <- 1 # remember: throw away draw
    #  }
    #}
    
    #if(ind == 1){
    #  next # throw away draw
    #}
    
    
    #########################################################################
    ######################### Check elasticity assumptions ##################
    #########################################################################
    
    ###### ELASTICITY RESTRICTION 1:  Bound price elasticity of gas supply after
    ######                            gas specific demand shock to 0.065
    ###### (This value represents the point estimate of Mason and Roberts, 2018,
    ######  for the 25% most productive active gas wells in Wyoming. The less
    ######  productive gas wells had way lower estimate, the estimate for all wells
    ######  being below 0.03. Our bound is therefore somehow prudent.)
    
    if((c1 %*% q3)/(c3 %*% q3) > ela){
      next # throw away draw
    }
    
    ###### ELASTICITY RESTRICTION 2:  Bound price elasticity of gas supply after
    ######                            aggregate demand shock to 0.065
    ###### Although this particular particular is somehow harder to motivate,
    ###### it is needed to rule out unplausible IRFs.
    
    if((c1 %*% q2)/(c3 %*% q2) > ela){ 
      next # throw away draw
    }
    
    ###### ELASTICITY RESTRICTION 3:  An one standard deviation gas demand shock
    ######                            is allowed to decrease industrial production
    ######                            by no more than 0.004 (0.008 is standard dev.
    ######                            of reduced form error for ipd)
    
    if(c2 %*% q3 < -a34){
      next # throw away draw
    }
    
    #########################################################################
    ################ If draw has been kept: compute IRFs and get their pdf ##
    #########################################################################
    
    # Fill result list
    resu$BET <- append(resu$BET, list(BET))
    resu$C <- append(resu$C, list(C))
    resu$Q <- append(resu$Q, list(Q))
    
    # empty vectors for IRFs
    rig_rig <- rep(NA, h+1)
    rig_gpd <- rep(NA, h+1)
    rig_ipd <- rep(NA, h+1)
    rig_rpg <- rep(NA, h+1)
    gpd_rig <- rep(NA, h+1)
    gpd_gpd <- rep(NA, h+1)
    gpd_ipd <- rep(NA, h+1)
    gpd_rpg <- rep(NA, h+1)
    ipd_rig <- rep(NA, h+1)
    ipd_gpd <- rep(NA, h+1)
    ipd_ipd <- rep(NA, h+1)
    ipd_rpg <- rep(NA, h+1)
    rpg_rig <- rep(NA, h+1)
    rpg_gpd <- rep(NA, h+1)
    rpg_ipd <- rep(NA, h+1)
    rpg_rpg <- rep(NA, h+1)
    
    for (j in 0:h) {
      tmp <- RIRF[[j+1]] %*% Q # structural IRF
      
      # add to vectors
      rig_rig[j+1] <- tmp[1,1]
      rig_gpd[j+1] <- tmp[2,1]
      rig_ipd[j+1] <- tmp[3,1]
      rig_rpg[j+1] <- tmp[4,1]
      gpd_rig[j+1] <- tmp[1,2]
      gpd_gpd[j+1] <- tmp[2,2]
      gpd_ipd[j+1] <- tmp[3,2]
      gpd_rpg[j+1] <- tmp[4,2]
      ipd_rig[j+1] <- tmp[1,3]
      ipd_gpd[j+1] <- tmp[2,3]
      ipd_ipd[j+1] <- tmp[3,3]
      ipd_rpg[j+1] <- tmp[4,3]
      rpg_rig[j+1] <- tmp[1,4]
      rpg_gpd[j+1] <- tmp[2,4]
      rpg_ipd[j+1] <- tmp[3,4]
      rpg_rpg[j+1] <- tmp[4,4]
    }
    
    #rig_gpd <- cumsum(rig_gpd)
    #gpd_gpd <- cumsum(gpd_gpd)
    #ipd_gpd <- cumsum(ipd_gpd)
    #rpg_gpd <- cumsum(rpg_gpd)
    
    # add to data frames
    resu$irig$rrig <- data.frame(resu$irig$rrig, rig_rig)
    resu$irig$rgpd <- data.frame(resu$irig$rgpd, rig_gpd)
    resu$irig$ripd <- data.frame(resu$irig$ripd, rig_ipd)
    resu$irig$rrpg <- data.frame(resu$irig$rrpg, rig_rpg)
    resu$igpd$rrig <- data.frame(resu$igpd$rrig, gpd_rig)
    resu$igpd$rgpd <- data.frame(resu$igpd$rgpd, gpd_gpd)
    resu$igpd$ripd <- data.frame(resu$igpd$ripd, gpd_ipd)
    resu$igpd$rrpg <- data.frame(resu$igpd$rrpg, gpd_rpg)
    resu$iipd$rrig <- data.frame(resu$iipd$rrig, ipd_rig)
    resu$iipd$rgpd <- data.frame(resu$iipd$rgpd, ipd_gpd)
    resu$iipd$ripd <- data.frame(resu$iipd$ripd, ipd_ipd)
    resu$iipd$rrpg <- data.frame(resu$iipd$rrpg, ipd_rpg)
    resu$irpg$rrig <- data.frame(resu$irpg$rrig, rpg_rig)
    resu$irpg$rgpd <- data.frame(resu$irpg$rgpd, rpg_gpd)
    resu$irpg$ripd <- data.frame(resu$irpg$ripd, rpg_ipd)
    resu$irpg$rrpg <- data.frame(resu$irpg$rrpg, rpg_rpg)
    
    rm(j, rig_rig, rig_gpd, rig_ipd, rig_rpg, gpd_rig, gpd_gpd, gpd_ipd, gpd_rpg, 
       ipd_rig, ipd_gpd, ipd_ipd, ipd_rpg, rpg_rig, rpg_gpd, rpg_ipd, rpg_rpg)
    
    # Dimension of VAR
    n <- 4
    
    # get posterior pdf of irf
    f <- irfpdf_blockrec(Q, n, p, BET, SIG, betbarT, VT, nuT, ST)
    resu$poster_f <- c(resu$poster_f, f)
    
    rm(f, n)
  }
  
  return(resu)
}