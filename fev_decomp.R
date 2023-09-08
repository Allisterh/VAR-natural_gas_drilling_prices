##############################################################################
### Script: FEV DECOMP
### Author: Valentin Winkler
### Date:   08/2022
### Description:  Compute forecast error variance decomposition for structural 
###               VAR model.
##############################################################################
### Inputs:   A     n x pn Matrix with reduced form slope parameters
###           C     Cholesky decomposition of estimated variance matrix
###           Q     drawn orthogonal matrix for sign identified models
###                 (default identity matrix)
###           h_fev Horizon of FEVD
###           plot  logical: should plot of FEVD be returned?
###           
### Outputs:  n Data Frames with FEVD and n plots
##############################################################################

# packages
library(ggplot2)
library(dplyr)

fev_decomp <- function(A, C, Q, h_fev, plot = T){
  # Get some parameters
  n <- dim(A)[1] # number of endogenous variables
  p <- dim(A)[2] / n # number of VAR-lags
  
  ############ Compute h-1 IRFs
  
  # construct companion matrix for VAR model
  AA <- matrix(0, n*p, n*p)
  AA[1:n,] <- A
  AA[-c(1:n),1:((p-1)*n)] <- diag((p-1)*n)
  
  # J-Matrix to get IRF as J (AA)^i J' C Q
  J <- matrix(0,n,n*p)
  J[1:n,1:n] <- diag(n)
  
  RIRF <- list(C %*% Q, J %*% AA %*% t(J) %*% C %*% Q)
  tmp <- AA
  
  for (i in 2:(h_fev-1)) {
    tmp <- tmp %*% AA
    RIRF <- append(RIRF, list(J %*% tmp %*% t(J) %*% C %*% Q))
  }
  
  ## result list
  resu <- list()
  nam <- c()
  
  for (i in 1:n) {
    resu_tmp <- data.frame(Horizon = 1:h_fev)
    name_tmp <- c("Horizon")
    for (j in 1:n) {
      err <- rep(0,h_fev)
      
      # update name vector
      name_tmp <- c(name_tmp, paste("contrib_var",j,sep = ""))
      
      for (v in 1:h_fev) {
        err[v] <- (RIRF[[v]][i,j])^2
      }
      # cumulate squared contributions
      err <- cumsum(err)
      resu_tmp <- data.frame(resu_tmp,  err)
    }
    
    # name data frame
    names(resu_tmp) <- name_tmp
    
    # normalize sum of forecast error contributions to 100
    for (v in 1:h_fev) {
      resu_tmp[v,-1] <- 100* resu_tmp[v,-1]/sum(resu_tmp[v,-1])
    }
    
    # Update result list
    resu <- append(resu, list(resu_tmp))
    
    # Update name
    nam <- c(nam, paste("decomp_var",i,sep=""))
  }
  
  names(resu) <- nam
  
  ######################## Plot FEVD for all variables #####################
  plots <- list()
  nam_plot <- c()
  
  for (i in 1:n) {
    D <- resu[[i]]
    
    pl <- pivot_longer(D, !"Horizon") %>%
      ggplot(aes(x = Horizon, y = value, fill = name),fig(5,12)) +
      geom_area() +
      xlab('Months') +
      ylab('Percent') +
      theme_bw()
    
    plots <- append(plots, list(pl))
    nam_plot <- c(nam_plot, paste("Plot_var",i,sep=""))
  }
  
  names(plots) <- nam_plot
  resu <- append(resu,plots)
  
  return(resu)
}

