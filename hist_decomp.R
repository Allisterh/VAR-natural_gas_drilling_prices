################################################################################
### Script: Historical Decomposition
### Author: Valentin Winkler
### Date:   08/2022
### Description:  Compute the historical decomposition of an n-dimensional VAR model.
###               Note that this scripts departs from Kilian and Lütkepohl (2018)
###               in that we do not rely on the MA(infinity) representation of
###               stable AR(p)-systems.
###               In particular we report the value y_t^i, which is y_t in the
###               case of all structural shocks apart from shock i being zero from t=1 on
################################################################################

######## Input
# Y        (n x T) Vector of endogenous variables
# Z        (M x T) Regressor matrix with columns [1, y_{t-1}',...,y_{t-p}',x_t']'
# B        (n x M) Slope coefficient matrix of the form [d B_1 ... B_p B_x] (reduced form)
# A_0^{-1} (n x n) inverse rotation matrix
# p        lag order

######## Output
# n (n x T) matrices Y^i, i = 1,..,n of historical decompositions
# 1 (n x T) matrix Y^c for comparison

##############################################################################

### load packages
library(dplyr)
library(tidyr)

hist_decomp <- function(Y, Z, B, A0_1, p){
  # get n and T
  n <- dim(Y)[1]
  T <- dim(Y)[2]
  M <- dim(B)[2]
  
  # get reduced form shocks
  E <- Y - B %*% Z
  
  # get A_0 matrix
  A_0 <- solve(A0_1)
  
  # get structural shocks
  EPS <- A_0 %*% E
  
  # result list and list names
  res <- list()
  nam <- c()
  
  # get decomposition for shock i
  for (i in 1:n) {
    
    # contrafactual structural and reduced form shocks
    EPS_i <- EPS*0
    EPS_i[i,] <- EPS[i,]
    E_i <- A0_1 %*% EPS_i
    
    # initialize Y_i and Z_i
    Y_i <- Y * 0
    Z_i <- Z
    
    # compute Y_t^i for each t
    for (t in 1:T) {
      Y_i[ ,t] <- B %*% Z_i[ ,t] + E_i[,t]
      
      # insert Y_t^i instead of Y_t in regressor matrix
      for (j in 1:p) {
        if(t + j <= T){
          Z_i[1 + (j-1)*n + (1:n),t+j] <- Y_i[ ,t]
        }
      }
    }
    
    res <- append(res, list(Y_i))
    nam <- c(nam, paste("Shock",i,sep = "_"))
  }
  
  ## get comparison values for each shock set to zero
  
  # zero shocks and initialize Y^c as well as Z^c
  E_c <- E * 0
  Y_c <- Y
  Z_c <- Z
  
  # compute Y_t^c for each t
  for (t in 1:T) {
    Y_c[ ,t] <- B %*% Z_c[ ,t] + E_c[,t]
    
    # insert Y_t^c instead of Y_t in regressor matrix
    for (j in 1:p) {
      if(t + j <= T){
        Z_c[1 + (j-1)*n + (1:n),t+j] <- Y_c[ ,t]
      }
    }
  }
  
  res <- append(res, list(Y_c))
  nam <- c(nam, "Comparison")
  
  names(res) <- nam
  return(res)
}


