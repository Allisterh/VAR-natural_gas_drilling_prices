################## Multivariate Gamma Function ########################
##### compute in log levels to avoid too big numbers

mult_gamma_log <- function(a,p){
  tmp <- log(pi^(p*(p-1)/4))
  
  for(i in 1:p){
    tmp <- tmp + log(gamma(a + (1-i)/2))
  }
  
  return(tmp)
}