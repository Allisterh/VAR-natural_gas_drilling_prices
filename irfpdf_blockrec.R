##############################################################################
### Script: IRF PDF
### Author: Valentin Winkler 
### Date:   08/2022
### Description:  Compute posterior value of IRF for a Bayesian SVAR model
###               with zero restrictions in the first row of A_0^-1 such that
###               the first structural shock is point identified and all else
###               set identified.
###               The model has to be estimated using sub-rotations and all
###               other model specifications like in Inoue and Kilian (2013).
##############################################################################
### Inputs:   Q       Orthogonal matrix satisfying sign restrictions
###           n       Dimensions of VAR-Model
###           p       lag order of VAR-Model
###           BET     drawn reduced form parameters
###           SIG     drawn VC-Matrix for errors
###           betbarT posterior mean of vec(BET)
###           VT      V1 + Z*Z'
###           nuT     df for inverse-Wishart distribution
###           ST      Scale matrix for inverse-Wishart distribution
###
### Outputs:  f       posterior value of IRF (in logs to avoid big/small number problems)
##############################################################################

source("mult_gamma_log.R") # multivariate gamma function
source("matrix_calc.R") # usefull functions for matrix calculus

irfpdf_blockrec <- function(Q, n, p, BET, SIG, betbarT, VT, nuT, ST){
  # compute some variables
  C <- t(chol(SIG))
  B0 <- C %*% Q
  Qstar <- Q[-1,-1]
  
  # Make sure that det(Qstar) = 1
  if(det(Qstar) < 0){
    W <- diag(c(-1,rep(1,n-2)))
    Qstar <- W %*% Qstar
    
    W2 <- diag(c(1,-1,rep(1,n-2)))
    Q <- W2 %*% Q
    C <- C %*% W2
    # B0 = C * W2 * W2 * Q remains unchanged
  }
  
  # Compute matrix K
  K <- matrix(0, nrow = n, ncol = n-1)
  K[-1,] <- diag(n-1)
  
  ########### (n-1)^2 x (n-1)(n-2)/2 Duplication matrix Dbar_1 that satisfies vec(S) = Dbar_1 s 
  Dbar_1 <- Dbarn(n-1)
  
  ########### Matrix I_{n^2} - L_n'L_n as needed in Kilian and Inoue (2018)
  ILL <- I_LL(n)
  
  
  # We now rely on a formula adapted from the corrigendum to Kilian and Inoue (2013)
  # published in 2018.
  # There it is proven that f = C1 * C2 * C3 * C4, where C1 is the determinant of a jacobian
  # for change of variables, C2 is the cond. density of slope parameters B, C3 is the density
  # of the VC-matrix and C4 is the density of s^*.
  # We report the density in log levels, i.e. we compute log(C1) + log(C2) + log(C3) + log(C4)
  
  ############### Compute C1 ################################################
  # This is the simplified inverse determinant of the jacobian associated with
  # the change of variables.
  
  C1 <- ((n^2 - n + 2)/2)*log(2) - n*p*log(abs(det(C)))
  
  for (i in 1:n) {
    C1 <- C1 + (n-i+1)*log(abs(C[i,i]))
  }
  
  C1 <- C1 -(n-2)*log(abs(det( diag(n-1) + Qstar  )))
  
  # complicated matrix for determinant in the middle of the formula
  tmpmat <- t(Dbar_1) %*% (t(K) %x% (t(K) %*% t(C))) %*% ILL %*% (K %x% (C %*% K)) %*% Dbar_1
  
  C1 <- C1 -(1/2)*log(abs(det(tmpmat)))
  
  ############### Compute C2 ################################################
  # C2 = f(A|SIG), conditional posterior density of the AR slope-parameters
  # Remember: Normal Distribution
  
  # get vec(A), A being estimated slope parameters
  A <- BET[,2:(p*n + 1)]
  a <- matrix(A, ncol = 1)
  
  # get VC-Matrix of vec(A)
  SIGBet <- solve(VT) %x% SIG
  SIGa <- SIGBet[(n+1):(p*n^2 + n), (n+1):(p*n^2 + n)]
  
  # get vec(Abar), Abar being expected value of slope parameter
  abar <- betbarT[(n+1):(p*n^2 + n) ,1]
  
  C2 <- ((2*pi)^n * det(SIGa))^(-1/2)*exp(-1/2 * (t(a - abar) %*% solve(SIGa) %*% (a - abar)))
  C2 <- log(C2)
  
  ############### Compute C3 ################################################
  # C3 = f(SIG), unconditional posterior density of innovation variance matrix
  # Remember: Inverse-Wishart Distribution
  
  # compute in log levels to avoid maximal values in R
  C3 <- (nuT/2)*log(det(ST)) - (nuT*p /2)*log(2) - mult_gamma_log(nuT/2, p) - ((nuT + p + 1)/2)*log(det(SIG)) - (1/2)*sum(diag(ST %*% solve(SIG)))
  
  # back to normal levels
  
  ############### Compute C4 ################################################
  # C4 = f(s^*), with s being the (n-1)(n-2)/2 below-diagonal elements of S^* and
  # S^* = I_{n-1} - 2(I_{n-1} + Qstar)^-1
  # the density of s^* can be computed according to a formula of León et al (2006)
  
  Sstar <- diag(n-1) - 2 * solve(diag(n-1) + Qstar)
  
  fac <- 1
  
  for (v in 2:(n-1)) {
    fac <- fac * gamma(v/(n-1))/(pi^(v/2))
  }
  
  C4 <- fac * 2^((n-2)*(n-3)/2) / det(diag(n-1) + Sstar)^(n-2)
  C4 <- log(C4)
  
  ### Compose and return posterior density
  f <- C1 + C2 + C3 + C4
  
  return(f)
}