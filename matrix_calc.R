##############################################################
######### Some useful functions for matrix calculus ##########
##############################################################
######### Author:   Valentin Winkler                ##########
######### Date:     08/2022                         ##########
##############################################################

#################### Special Matrices #################################

# Input:  dimension n
# Output: special n^2 x n(n-1)/2 Duplication matrix as in Kilian and Inoue (2018)
Dbarn <- function(n){
  res <- matrix(0, nrow = n^2, ncol = n*(n-1)/2)
  
  for (j in 1:(n-1)) {
    for (i in (j+1):n) {
      # define nxn matrices E_ij, E_ji
      Eij <- matrix(0, nrow = n, ncol = n)
      Eji <- matrix(0, nrow = n, ncol = n)
      
      Eij[i,j] <- 1
      Eji[j,i] <- 1
      
      # add to result
      res <- res + vec(Eij - Eji) %*% t(veck(Eij))
    }
  }
  return(res)
}

# Input:  dimension n
# Output: Matrix I_{n^2} - L_n'L_n as needed in Kilian and Inoue (2018)
I_LL <- function(n){
  res <- matrix(0 , nrow = n^2, ncol = n^2)
  
  for (j in 1:(n-1)) {
    for (i in (j+1):n) {
      # define nxn matrices E_ij, E_ji
      Eii <- matrix(0, nrow = n, ncol = n)
      Ejj <- matrix(0, nrow = n, ncol = n)
      
      Eii[i,i] <- 1
      Ejj[j,j] <- 1
      
      # add to result
      res <- res + Eii %x% Ejj
    }
  }
  return(res)
}

# Input:  dimension n
# Ouptut: special n(n-1) x n^2 matrix Q that contains only zeros and ones and
#         satisfies RR' = I and R'R = I - L_n'L_n as defined above

Rmat <- function(n){
  R <- matrix(nrow = n*(n-1), ncol = n^2)
    
  for (j in 1:(n-1)) {
    for (i in (j+1):n) {
      # define nxn matrices E_ij, E_ji
      Eij <- matrix(0, nrow = n, ncol = n)
      Eji <- matrix(0, nrow = n, ncol = n)
      
      Eij[i,j] <- 1
      Eji[j,i] <- 1
      
      R <- R + veck(Eij) %*% t(vec(Eji))
    }
  }
}

#################### Matrix to Vector #################################

###### Function VEC
# Input: nxm matrix M
# Output: nm x 1 vec(M) ... stacked columns of M
vec <- function(M){
  return(matrix(M, ncol = 1))
}

###### Function VECH
# Input: nxn matrix M
# Output: (n+1)n/2 x 1 vector vech(M) with below- and on-diagonal elements 
#         stacked column-wise
vech <- function(M){
  return(matrix(M[lower.tri(M, diag = T)], ncol = 1))
}

###### Function VECK
# Input: nxn matrix M
# Output: (n-1)n/2 x 1 vector veck(M) with below- and NOT on-diagonal elements 
#         stacked column-wise
veck <- function(M){
  return(matrix(M[lower.tri(M, diag = F)], ncol = 1))
}
