##############################################################################
### Script: Bayes 4x4 VAR
### Author: Valentin Winkler
### Date:   08/2023
### Description:  Functions that take draws of structural impulse response functions
###               as provided by 'bayes_4_4_var.R' and return posterior estimates
###               of price elasticity of drilling and demand side contribution to gas
###               price collapse, respectively.


######### Function drill_rob #######################################
f_elast_rob <- function(mod, al = 0.32){
  ### Input:   model object from 'bayes_4_4_var.R' and confidence level 'al'
  ### Ouptput: lower bound, median and upper bound (corresponding to al) for posterior distributions
  
  res <- rep(0,32)
  
  # difference reaction of rig to ipd and gpd shock
  r_ipd <- as.matrix(mod$iipd$rrig[,-1], nrow = 13)
  p_ipd <-  matrix(1, nrow = 13, ncol = 1) %*%  as.matrix(mod$iipd$rrpg[1,-1], nrow = 1)
  r_ipd <- r_ipd/p_ipd
  
  r_gpd <- as.matrix(mod$igpd$rrig[,-1], nrow = 13)
  p_gpd <- matrix(1, nrow = 13, ncol = 1) %*% as.matrix(mod$igpd$rrpg[1,-1], nrow = 1)
  r_gpd <- r_gpd/p_gpd
  
  diff_ipd_gpd <- r_ipd - r_gpd
  
  diff_ipd_gpd2 <- t(diff_ipd_gpd[c(4,7,10,13),]) %>%
    as.data.frame()
  
  names(diff_ipd_gpd2) <- c("h3", "h6", "h9", "h12")
  
  # difference reaction of rig to ipd and rpg shock
  r_rpg <- as.matrix(mod$irpg$rrig[,-1], nrow = 13)
  p_rpg <- matrix(1, nrow = 13, ncol = 1) %*% as.matrix(mod$irpg$rrpg[1,-1], nrow = 1)
  r_rpg <- r_rpg/p_rpg
  
  diff_ipd_rpg <- r_ipd - r_rpg
  
  diff_ipd_rpg2 <- t(diff_ipd_rpg[c(4,7,10,13),]) %>%
    as.data.frame()
  
  names(diff_ipd_rpg2) <- c("h3", "h6", "h9", "h12")
  
  # fill results vector
  res[1:3]   <- quantile(diff_ipd_gpd2$h3, probs = c(al/2, 0.5, 1-al/2))
  res[4]     <- round(mean(diff_ipd_gpd2$h3 < 0), digits = 3)
  res[5:7]   <- quantile(diff_ipd_gpd2$h6, probs = c(al/2, 0.5, 1-al/2))
  res[8]     <- round(mean(diff_ipd_gpd2$h6 < 0), digits = 3)
  res[9:11]  <- quantile(diff_ipd_gpd2$h9, probs = c(al/2, 0.5, 1-al/2))
  res[12]    <- round(mean(diff_ipd_gpd2$h9 < 0), digits = 3)
  res[13:15] <- quantile(diff_ipd_gpd2$h12, probs = c(al/2, 0.5, 1-al/2))
  res[16]    <- round(mean(diff_ipd_gpd2$h12 < 0), digits = 3)
  res[17:19] <- quantile(diff_ipd_rpg2$h3, probs = c(al/2, 0.5, 1-al/2))
  res[20]    <- round(mean(diff_ipd_rpg2$h3 < 0), digits = 3)
  res[21:23] <- quantile(diff_ipd_rpg2$h6, probs = c(al/2, 0.5, 1-al/2))
  res[24]    <- round(mean(diff_ipd_rpg2$h6 < 0), digits = 3)
  res[25:27] <- quantile(diff_ipd_rpg2$h9, probs = c(al/2, 0.5, 1-al/2))
  res[28]    <- round(mean(diff_ipd_rpg2$h9 < 0), digits = 3)
  res[29:31] <- quantile(diff_ipd_rpg2$h12, probs = c(al/2, 0.5, 1-al/2))
  res[32]    <- round(mean(diff_ipd_rpg2$h12 < 0), digits = 3)
  
  return(res)
}

######### Function contrib_rob #######################################
f_contrib_rob <- function(mod, p, s1, e1, s2, e2, s3, e3, al = 0.32){
  ### Input:   model object from 'bayes_4_4_var.R', start and end date of three time windows and confidence level 'al'
  ### Ouptput: lower bound, median and upper bound (corresponding to al) for posterior distributions
  
  res <- rep(0,9)
  
  # get variables
  Y <- mod$LS_model$Y
  Z <- mod$LS_model$Z
  r <- length(mod$poster_f)
  all <- Y[4,]
  
  # get gas price decrease in 3 windows
  l1 <- all[e1] - all[s1]
  l2 <- all[e2] - all[s2]
  l3 <- all[e3] - all[s3]
  
  # initialize result vecotrs
  win1 <- mod$poster_f * 0
  win2 <- mod$poster_f * 0
  win3 <- mod$poster_f * 0
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = r, initial = 0, style = 3)
  
  # get fractions for every model
  for (i in 1:r) {
    # update progress bar
    setTxtProgressBar(pb,i)
    
    B <- mod$BET[[i]]
    A0_1 <- mod$C[[i]] %*% mod$Q[[i]]
    
    hi <- hist_decomp(Y, Z, B, A0_1, p)
    ipd <- hi$Shock_3[4,]
    rpg <- hi$Shock_4[4,]
    com <- hi$Comparison[4,]
    
    # episode 1
    ea <- (ipd[e1] - com[e1]) - (ipd[s1] - com[s1])
    de <- (rpg[e1] - com[e1]) - (rpg[s1] - com[s1])
    win1[i] <- ((ea + de)/l1) * 100
    
    # episode 2
    ea <- (ipd[e2] - com[e2]) - (ipd[s2] - com[s2])
    de <- (rpg[e2] - com[e2]) - (rpg[s2] - com[s2])
    win2[i] <- ((ea + de)/l2) * 100
    
    # episode 3
    ea <- (ipd[e3] - com[e3]) - (ipd[s3] - com[s3])
    de <- (rpg[e3] - com[e3]) - (rpg[s3] - com[s3])
    win3[i] <- ((ea + de)/l3) * 100
  }
  
  # fill results vector
  res[1:3]   <- quantile(win1, probs = c(al/2, 0.5, 1-al/2))
  res[4:6]   <- quantile(win2, probs = c(al/2, 0.5, 1-al/2))
  res[7:9]   <- quantile(win3, probs = c(al/2, 0.5, 1-al/2))
  
  close(pb) # close progress bar
  return(res)
}



