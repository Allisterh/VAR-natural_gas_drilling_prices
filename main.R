##############################################################################
### Script: 4x4 Gasmodels
### Author: Valentin Winkler
### Date:   08/2022
### Description:  Estimate sign and exclusion restricted, bayesian SVAR-models 
###               for the US gas market
##############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(vars)
library(ggpubr)
library(reshape2)
library(extrafont)
library(scales)
source("bayes_4_4_var.R")
source("fev_decomp.R")
source("hist_decomp.R")
source("robustness_functions.R")

loadfonts()

# for (1-alpha) highest density band
alph = 0.32

# horizon of fevd
h_fev <- 14

# Draws of Reduced Form Parameters
M <- 500
# Draws of Orthogonal Matrix for each reduced form draw
K <- 10000
p <- 6

set.seed(11811850)

###############################################################################################
############################## Estimate Model and get IRFs ####################################
###############################################################################################

# Estimate VAR for full sample until before covid-shock
full_sample <- bayes_4_4_var(star = as.Date("1993-11-01"), en = as.Date("2020-01-01"), M, K)
#save(full_sample, file = "moddata/fulldata.RData")

# Sort draws according to posterior values

ord <- order(full_sample$poster_f, decreasing = T)

####################### Impact rig
rig_rig <- full_sample$irig$rrig
rig_rig[,-1] <- rig_rig[ ,1 + ord]
na <- names(rig_rig)
na[2] <- "modal"
names(rig_rig) <- na

rig_gpd <- full_sample$irig$rgpd
rig_gpd[,-1] <- rig_gpd[ ,1 + ord]
na <- names(rig_gpd)
na[2] <- "modal"
names(rig_gpd) <- na

rig_ipd <- full_sample$irig$ripd
rig_ipd[,-1] <- rig_ipd[ ,1 + ord]
na <- names(rig_ipd)
na[2] <- "modal"
names(rig_ipd) <- na

rig_rpg <- full_sample$irig$rrpg
rig_rpg[,-1] <- rig_rpg[ ,1 + ord]
na <- names(rig_rpg)
na[2] <- "modal"
names(rig_rpg) <- na

####################### Impact gpd
gpd_rig <- full_sample$igpd$rrig
gpd_rig[,-1] <- gpd_rig[ ,1 + ord]
na <- names(gpd_rig)
na[2] <- "modal"
names(gpd_rig) <- na

gpd_gpd <- full_sample$igpd$rgpd
gpd_gpd[,-1] <- gpd_gpd[ ,1 + ord]
na <- names(gpd_gpd)
na[2] <- "modal"
names(gpd_gpd) <- na

gpd_ipd <- full_sample$igpd$ripd
gpd_ipd[,-1] <- gpd_ipd[ ,1 + ord]
na <- names(gpd_ipd)
na[2] <- "modal"
names(gpd_ipd) <- na

gpd_rpg <- full_sample$igpd$rrpg
gpd_rpg[,-1] <- gpd_rpg[ ,1 + ord]
na <- names(gpd_rpg)
na[2] <- "modal"
names(gpd_rpg) <- na

####################### Impact ipd
ipd_rig <- full_sample$iipd$rrig
ipd_rig[,-1] <- ipd_rig[ ,1 + ord]
na <- names(ipd_rig)
na[2] <- "modal"
names(ipd_rig) <- na

ipd_gpd <- full_sample$iipd$rgpd
ipd_gpd[,-1] <- ipd_gpd[ ,1 + ord]
na <- names(ipd_gpd)
na[2] <- "modal"
names(ipd_gpd) <- na

ipd_ipd <- full_sample$iipd$ripd
ipd_ipd[,-1] <- ipd_ipd[ ,1 + ord]
na <- names(ipd_ipd)
na[2] <- "modal"
names(ipd_ipd) <- na

ipd_rpg <- full_sample$iipd$rrpg
ipd_rpg[,-1] <- ipd_rpg[ ,1 + ord]
na <- names(ipd_rpg)
na[2] <- "modal"
names(ipd_rpg) <- na

####################### Impact rpg
rpg_rig <- full_sample$irpg$rrig
rpg_rig[,-1] <- rpg_rig[ ,1 + ord]
na <- names(rpg_rig)
na[2] <- "modal"
names(rpg_rig) <- na

rpg_gpd <- full_sample$irpg$rgpd
rpg_gpd[,-1] <- rpg_gpd[ ,1 + ord]
na <- names(rpg_gpd)
na[2] <- "modal"
names(rpg_gpd) <- na

rpg_ipd <- full_sample$irpg$ripd
rpg_ipd[,-1] <- rpg_ipd[ ,1 + ord]
na <- names(rpg_ipd)
na[2] <- "modal"
names(rpg_ipd) <- na

rpg_rpg <- full_sample$irpg$rrpg
rpg_rpg[,-1] <- rpg_rpg[ ,1 + ord]
na <- names(rpg_rpg)
na[2] <- "modal"
names(rpg_rpg) <- na

################################### Plot IRFs with 68% highest density band
obs <- dim(rig_rig)[2] - 1
obs_band <- ceiling((1 - alph)*obs)

options(digits = 4)

########## Response to drilling shock
p11 <- rig_rig[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "forest green") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  geom_hline(aes(yintercept = 0), col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Rig Count') +
  ggtitle("Drilling Shock") +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ylim(-0.04,0.1)+ 
  theme(text= element_text(family="LM Roman 10")) 

p21 <- rig_gpd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "forest green") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  geom_hline(aes(yintercept = 0), col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Gas Production') +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ylim(-0.02,0.011)+ 
  theme(text= element_text(family="LM Roman 10"))

p31 <- rig_ipd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "forest green") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Industrial Production')+
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.018,0.018)+ 
  theme(text= element_text(family="LM Roman 10"))

p41 <- rig_rpg[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "forest green") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('Months after Shock') +
  ylab('Real Price of Gas')+
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.06, 0.13)+ 
  theme(text= element_text(family="LM Roman 10"))

########## Response to gas supply shock
p12 <- gpd_rig[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "cornflowerblue") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  geom_hline(aes(yintercept = 0), col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Rig Count') +
  ggtitle("Gas Supply Shock") +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.04,0.1)+ 
  theme(text= element_text(family="LM Roman 10"))

p22 <- gpd_gpd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "cornflowerblue") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  geom_hline(aes(yintercept = 0), col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Gas Production') +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ylim(-0.02,0.011)+ 
  theme(text= element_text(family="LM Roman 10"))

p32 <- gpd_ipd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "cornflowerblue") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Industrial Production')+
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.018,0.018)+ 
  theme(text= element_text(family="LM Roman 10"))

p42 <- gpd_rpg[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "cornflowerblue") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('Months after Shock') +
  ylab('Real Price of Gas')+
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.06, 0.13)+ 
  theme(text= element_text(family="LM Roman 10"))

######### Response to economic activity shock
p13 <- ipd_rig[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "violet") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  geom_hline(aes(yintercept = 0), col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Rig Count') +
  ggtitle("Economic Activity Shock") +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.04,0.1)+ 
  theme(text= element_text(family="LM Roman 10"))

p23 <- ipd_gpd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "violet") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Gas Production') +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.02,0.011)+ 
  theme(text= element_text(family="LM Roman 10"))

p33 <- ipd_ipd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "violet") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Industrial Production') +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.018,0.018)+ 
  theme(text= element_text(family="LM Roman 10"))

p43 <- ipd_rpg[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "violet") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('Months after Shock') +
  ylab('Real Price of Gas') +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.06, 0.13)+ 
  theme(text= element_text(family="LM Roman 10"))

########## Response to gas demand shock
p14 <- gpd_rig[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "red") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  geom_hline(aes(yintercept = 0), col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Rig Count') +
  ggtitle("Gas Demand Shock") +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.04,0.1)+ 
  theme(text= element_text(family="LM Roman 10"))

p24 <- rpg_gpd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "red") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Gas Production') +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.02,0.011)+ 
  theme(text= element_text(family="LM Roman 10"))

p34 <- rpg_ipd[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "red") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('') +
  ylab('Industrial Production') +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.018,0.018)+ 
  theme(text= element_text(family="LM Roman 10"))

p44 <- rpg_rpg[ ,1:(1+obs_band)] %>%
  melt(id.vars = c("Horizon", "modal")) %>%
  ggplot(aes(x = Horizon, y = value, group = variable)) +
  geom_line(size = 0.2, alpha = 1, col = "red") +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_line(aes(y = modal), size = 0.75, alpha = 1, col = "black") +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab('Months after Shock') +
  ylab('Real Price of Gas') +
  scale_x_continuous(breaks = pretty_breaks())+
  ylim(-0.06, 0.13)+ 
  theme(text= element_text(family="LM Roman 10"))

plot <- ggarrange(p11, p12, p13, p14, p21, p22, p23, p24, p31, p32, p33, p34, p41, p42, p43, p44, ncol = 4, nrow = 4, widths = 8, heights = 8)
plot
ggsave("modfig/test_44IRF_full_sample_bayes.png", device = "png", width = 10, height = 8)

#####################################################################################################
################################ Forecast Error Variance Decomposition ##############################
#####################################################################################################

# Find modal model
mod <- order(full_sample$poster_f, decreasing = T)[1]

# fevd model
A <- full_sample$BET[[mod]]
A <- A[,2:(4*p +1)]
C <- full_sample$C[[mod]]
Q <- full_sample$Q[[mod]]

fev_full <- fev_decomp(A, C, Q, h_fev, plot = T)

############## Format plots
colorinf <- scale_fill_manual(
  name = "Contributions",
  labels = c("Drilling", "Gas Supply", "Economic Activity", "Gas Demand"),
  values = c("forest green", "cornflowerblue", "grey", "red")
)

marg <- theme(plot.margin = unit(c(0,0.2,0,0.2),"cm")) +
  theme(text= element_text(family="LM Roman 10"))

p1 <- fev_full$Plot_var1 +
  ggtitle("Rig Count") +
  colorinf +
  xlab('') + marg

p2 <- fev_full$Plot_var2 +
  ggtitle("Gas Production") +
  colorinf +
  xlab('') + marg

p3 <- fev_full$Plot_var3 +
  ggtitle("Industrial Production") +
  colorinf +
  xlab('') + marg

p4 <- fev_full$Plot_var4 +
  ggtitle("Gas Price") +
  colorinf + marg

dev.new(width = 15, height = 20)
plot <- ggarrange(p1, p2, p3, p4, ncol = 1, common.legend = T, legend = "right")
plot
ggsave("modfig/44full_fevd.png", device = "png", width = 9, height = 5.5)

###################################################################################################################
######################################## Historical Decomposition #################################################
###################################################################################################################

### Full Sample
Y <- full_sample$LS_model$Y
Z <- full_sample$LS_model$Z

# order for modal model
mod <- order(full_sample$poster_f, decreasing = T)[1]

B <- full_sample$BET[[mod]]
A0_1 <- full_sample$C[[mod]] %*% full_sample$Q[[mod]]

hist_full <- hist_decomp(Y, Z, B, A0_1, p)

### Decomposition gas price

gp <- gp %>%
  mutate(all_d = all - Comparison,
         rig_d = rig - Comparison,
         gpd_d = gpd - Comparison,
         ipd_d = ipd - Comparison,
         rpg_d = rpg - Comparison,
         com_d = Comparison - Comparison)

nav <- rep(NA, dim(gp)[1])

pl_gp <- data.frame(Date = rep(gp$Date, 5),
                    Nam = rep(c("Real Price", "Drilling", "Supply", "Econ. Activity", "Gas Demand"),rep(dim(gp)[1],5)),
                    all = c(gp$all, gp$all_d, gp$all_d, gp$all_d, gp$all_d),
                    rig = c(nav, gp$rig_d, nav, nav, nav),
                    gpd = c(nav, nav, gp$gpd_d, nav, nav),
                    ipd = c(nav, nav, nav, gp$ipd_d, nav),
                    rpg = c(nav, nav, nav, nav, gp$rpg_d),
                    Comparison = c(gp$Comparison, gp$com_d, gp$com_d, gp$com_d, gp$com_d))

pl_gp_l <- pivot_longer(pl_gp, !c("Date", "Nam"))

pl_gp_l$Nam <- factor(pl_gp_l$Nam, levels = c("Real Price", "Drilling", "Supply", "Econ. Activity", "Gas Demand"))
pl_gp_l$name <- factor(pl_gp_l$name, levels = c("Comparison", "all", "rig", "gpd", "ipd", "rpg"))

ggplot(pl_gp_l, aes(x = Date, y = value, col = name, size = name)) +
  geom_line()+
  facet_grid(rows = vars(Nam), scales = "free_y") +
  scale_color_manual(values = c("black", "dark grey", "forest green", "cornflowerblue", "violet", "red")) +
  scale_size_manual(values = c(0.5, 0.7, 0.7, 0.7, 0.7, 0.7))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 6.5, color = "white", face = "bold"),
        strip.background = element_rect(color = "dark grey", fill = "dark grey")) +
  xlab('') +
  ylab('')+ 
  theme(text= element_text(family="LM Roman 10"))
  
ggsave("modfig/rpg_hist2.png", device = "png", width = 6.5, height = 4.5)

### Decomposition rig count

rg <- rg %>%
  mutate(all_d = all - Comparison,
         rig_d = rig - Comparison,
         gpd_d = gpd - Comparison,
         ipd_d = ipd - Comparison,
         rpg_d = rpg - Comparison,
         com_d = Comparison - Comparison)

nav <- rep(NA, dim(rg)[1])

pl_rg <- data.frame(Date = rep(rg$Date, 5),
                    Nam = rep(c("Rig Count", "Drilling", "Supply", "Econ. Activity", "Gas Demand"),rep(dim(rg)[1],5)),
                    all = c(rg$all, rg$all_d, rg$all_d, rg$all_d, rg$all_d),
                    rig = c(nav, rg$rig_d, nav, nav, nav),
                    gpd = c(nav, nav, rg$gpd_d, nav, nav),
                    ipd = c(nav, nav, nav, rg$ipd_d, nav),
                    rpg = c(nav, nav, nav, nav, rg$rpg_d),
                    Comparison = c(rg$Comparison, rg$com_d, rg$com_d, rg$com_d, rg$com_d))

pl_rg_l <- pivot_longer(pl_rg, !c("Date", "Nam"))

pl_rg_l$Nam <- factor(pl_rg_l$Nam, levels = c("Rig Count", "Drilling", "Supply", "Econ. Activity", "Gas Demand"))
pl_rg_l$name <- factor(pl_rg_l$name, levels = c("Comparison", "all", "rig", "gpd", "ipd", "rpg"))

ggplot(pl_rg_l, aes(x = Date, y = value, col = name, size = name)) +
  geom_line()+
  facet_grid(rows = vars(Nam), scales = "free_y") +
  scale_color_manual(values = c("black", "dark grey", "forest green", "cornflowerblue", "violet", "red")) +
  scale_size_manual(values = c(0.5, 0.7, 0.7, 0.7, 0.7, 0.7))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 6.5, color = "white", face = "bold"),
        strip.background = element_rect(color = "dark grey", fill = "dark grey")) +
  xlab('') +
  ylab('')+ 
  theme(text= element_text(family="LM Roman 10"))

ggsave("modfig/rig_hist2.png", device = "png", width = 6.5, height = 4.5)

### Plot decomposition for gas production
pr <- data.frame(Date = seq.Date(as.Date("1994-05-01"), length.out = dim(Y)[2], by = "month"),
                 all = Y[2,],
                 rig = hist_full$Shock_1[2,],
                 gpd = hist_full$Shock_2[2,],
                 ipd = hist_full$Shock_3[2,],
                 rpg = hist_full$Shock_4[2,],
                 Comparison = hist_full$Comparison[2,])

gp_l <- pivot_longer(pr, !c("Date", "Comparison"))
gp_l$name <- factor(gp_l$name, levels = c("all", "rig", "gpd", "ipd", "rpg"),
                    labels = c("Gas Production", "Drilling", "Supply", "Econ. Activity", "Gas Demand"))

ggplot(gp_l, aes(x = Date)) +
  geom_line(aes(y = Comparison), col = "black") +
  geom_line(aes(y = value, col = name)) +
  facet_grid(rows = vars(name)) +
  scale_color_manual(values = c("dark grey", "forest green", "cornflowerblue", "violet", "red")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 6.5, color = "white", face = "bold"),
        strip.background = element_rect(color = "dark grey", fill = "dark grey")) +
  xlab('') +
  ylab('')+ 
  theme(text= element_text(family="LM Roman 10"))

ggsave("modfig/gpd_hist.png", device = "png", width = 6.5, height = 4.5)

### Plot decomposition for industrial production
ip <- data.frame(Date = seq.Date(as.Date("1994-05-01"), length.out = dim(Y)[2], by = "month"),
                 all = Y[3,],
                 rig = hist_full$Shock_1[3,],
                 gpd = hist_full$Shock_2[3,],
                 ipd = hist_full$Shock_3[3,],
                 rpg = hist_full$Shock_4[3,],
                 Comparison = hist_full$Comparison[3,])

ip_l <- pivot_longer(ip, !c("Date", "Comparison"))
ip_l$name <- factor(ip_l$name, levels = c("all", "rig", "gpd", "ipd", "rpg"),
                    labels = c("Industrial Production", "Drilling", "Supply", "Econ. Activity", "Gas Demand"))

ggplot(ip_l, aes(x = Date)) +
  geom_line(aes(y = Comparison), col = "black") +
  geom_line(aes(y = value, col = name)) +
  facet_grid(rows = vars(name)) +
  scale_color_manual(values = c("dark grey", "forest green", "cornflowerblue", "violet", "red")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 6.5, color = "white", face = "bold"),
        strip.background = element_rect(color = "dark grey", fill = "dark grey")) +
  xlab('')+
  ylab('')+ 
  theme(text= element_text(family="LM Roman 10"))

ggsave("modfig/ipd_hist.png", device = "png", width = 6.5, height = 4.5)

##################################################################################################################
#################### Relative Contributions of shocks to gas price collapse ######################################
##################################################################################################################

# Windows: 2008-06-01 to 2009-09-01, 2008-06-01 to 2012-04-01 and 2007-06-01 to 2013-06-01
contr <- data.frame(Window = c("2008M06-2009M09", "2008M06-2012M04", "2007M06-2013M06"))

# change in gas price
contr$Price <- c(subset(gp, Date == as.Date("2009-09-01"))$all - subset(gp, Date == as.Date("2008-06-01"))$all,
                      subset(gp, Date == as.Date("2012-04-01"))$all - subset(gp, Date == as.Date("2008-06-01"))$all,
                 subset(gp, Date == as.Date("2013-06-01"))$all - subset(gp, Date == as.Date("2007-06-01"))$all)

# structural component
contr$Structural <- c(subset(gp, Date == as.Date("2009-09-01"))$Comparison - subset(gp, Date == as.Date("2008-06-01"))$Comparison,
                      subset(gp, Date == as.Date("2012-04-01"))$Comparison - subset(gp, Date == as.Date("2008-06-01"))$Comparison,
                      subset(gp, Date == as.Date("2013-06-01"))$Comparison - subset(gp, Date == as.Date("2007-06-01"))$Comparison)

# drilling
contr$rig <- c((subset(gp, Date == as.Date("2009-09-01"))$rig - subset(gp, Date == as.Date("2009-09-01"))$Comparison) 
                - (subset(gp, Date == as.Date("2008-06-01"))$rig - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2012-04-01"))$rig - subset(gp, Date == as.Date("2012-04-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$rig - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2013-06-01"))$rig - subset(gp, Date == as.Date("2013-06-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2007-06-01"))$rig - subset(gp, Date == as.Date("2007-06-01"))$Comparison))

# supply
contr$gpd <- c((subset(gp, Date == as.Date("2009-09-01"))$gpd - subset(gp, Date == as.Date("2009-09-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$gpd - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2012-04-01"))$gpd - subset(gp, Date == as.Date("2012-04-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$gpd - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2013-06-01"))$gpd - subset(gp, Date == as.Date("2013-06-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2007-06-01"))$gpd - subset(gp, Date == as.Date("2007-06-01"))$Comparison))

# economic activity
contr$ipd <- c((subset(gp, Date == as.Date("2009-09-01"))$ipd - subset(gp, Date == as.Date("2009-09-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$ipd - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2012-04-01"))$ipd - subset(gp, Date == as.Date("2012-04-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$ipd - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2013-06-01"))$ipd - subset(gp, Date == as.Date("2013-06-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2007-06-01"))$ipd - subset(gp, Date == as.Date("2007-06-01"))$Comparison))

# gas demand
contr$rpg <- c((subset(gp, Date == as.Date("2009-09-01"))$rpg - subset(gp, Date == as.Date("2009-09-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$rpg - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2012-04-01"))$rpg - subset(gp, Date == as.Date("2012-04-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2008-06-01"))$rpg - subset(gp, Date == as.Date("2008-06-01"))$Comparison),
               (subset(gp, Date == as.Date("2013-06-01"))$rpg - subset(gp, Date == as.Date("2013-06-01"))$Comparison) 
               - (subset(gp, Date == as.Date("2007-06-01"))$rpg - subset(gp, Date == as.Date("2007-06-01"))$Comparison))

contr <- contr %>%
  pivot_longer(!c("Window", "Price"))
  
contr$name <- factor(contr$name, levels = c("Structural", "rig", "gpd", "ipd", "rpg"))
contr$Window <- factor(contr$Window, levels = c("2008M06-2009M09", "2008M06-2012M04", "2007M06-2013M06"))

ggplot(contr, aes(x = Window, y = value, fill = name)) +
  geom_histogram(stat = "identity", position = "stack") +
  geom_hline(aes(yintercept = 0))+
  geom_segment(aes(x = 0.54, xend = 1.46, y = contr[1,]$Price, yend = contr[1,]$Price), size = 1.8) +
  geom_segment(aes(x = 1.54, xend = 2.46, y = contr[6,]$Price, yend = contr[6,]$Price), size = 1.8) +
  geom_segment(aes(x = 2.54, xend = 3.46, y = contr[11,]$Price, yend = contr[11,]$Price), size = 1.8) +
  scale_fill_manual(name = "Contribution",
                    values = c("dark grey", "forest green", "cornflowerblue", "violet", "red"),
                    labels = c("Trend", "Drilling", "Gas Supply", "Econ. Activity", "Gas Demand")) +
  xlab('')+
  ylab('')+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 10))+ 
  theme(text= element_text(family="LM Roman 10"))

ggsave("modfig/contributions.png", device = "png", width = 6.5, height = 3)

# calculate fraction of demand driven effects for 2007M06-2013M06
dem <- abs(contr[15,4]+contr[14,4])
all <- dem + abs(contr[13,4]+contr[12,4]+contr[11,4])

##################################################################################################
################# Posterior distribution of fraction of gas price collapse #######################
######################### that is explained by demand side shocks ################################
##################################################################################################

r <- length(full_sample$poster_f)
D <- seq.Date(as.Date("1994-05-01"), length.out = dim(Y)[2], by = "month")
all <- Y[4,]

# start and end date index of windows
s1 <- which(D == as.Date("2008-06-01"))
e1 <- which(D == as.Date("2009-09-01"))
s2 <- which(D == as.Date("2008-06-01"))
e2 <- which(D == as.Date("2012-04-01"))
s3 <- which(D == as.Date("2007-06-01"))
e3 <- which(D == as.Date("2013-06-01"))

# get gas price decrease in 3 windows
l1 <- all[e1] - all[s1]
l2 <- all[e2] - all[s2]
l3 <- all[e3] - all[s3]

# initialize result vecotrs
win1 <- full_sample$poster_f * 0
win2 <- full_sample$poster_f * 0
win3 <- full_sample$poster_f * 0

# get fractions for modal models
mod_win1 <- 100*(subset(contr, Window == "2008M06-2009M09" & name == "ipd")$value + 
                subset(contr, Window == "2008M06-2009M09" & name == "rpg")$value)/l1
mod_win2 <- 100*(subset(contr, Window == "2008M06-2012M04" & name == "ipd")$value + 
               subset(contr, Window == "2008M06-2012M04" & name == "rpg")$value)/l2
mod_win3 <- 100*(subset(contr, Window == "2007M06-2013M06" & name == "ipd")$value + 
               subset(contr, Window == "2007M06-2013M06" & name == "rpg")$value)/l3

for (i in 1:r) {
  B <- full_sample$BET[[i]]
  A0_1 <- full_sample$C[[i]] %*% full_sample$Q[[i]]
  
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

pl1 <- data.frame(win = win1) %>%
ggplot(aes(x = win)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = quantile(win1, probs = 0.16)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(win1, probs = 0.84)), linetype = "dashed") +
  geom_vline(aes(xintercept = mod_win1), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("2008M06-2009M09")+
  xlab("Percent")+
  xlim(0,120)+ 
  theme(text= element_text(family="LM Roman 10"))

pl2 <- data.frame(win = win2) %>%
  ggplot(aes(x = win)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = mod_win2), col = "red")+
  geom_vline(aes(xintercept = quantile(win2, probs = 0.16)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(win2, probs = 0.84)), linetype = "dashed") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("2008M06-2012M04")+
  xlab("Percent")+
  ylab("")+
  xlim(0,120)+ 
  theme(text= element_text(family="LM Roman 10"))

pl3 <- data.frame(win = win3) %>%
  ggplot(aes(x = win)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = quantile(win3, probs = 0.16)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(win3, probs = 0.84)), linetype = "dashed") +
  geom_vline(aes(xintercept = mod_win3), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("2007M06-2013M06")+
  xlab("Percent")+
  ylab("")+
  xlim(0,120)+ 
  theme(text= element_text(family="LM Roman 10"))

ggarrange(pl1, pl2, pl3, ncol = 3)
ggsave("modfig/post_contr.png", device = "png", width = 7.5, height = 2.25)

##################################################################################################
#################### IRF elasticities rig reaction to ipd - gpd and ipd - rpg ####################
##################################################################################################

######### difference reaction of rig to ipd and gpd shock
r_ipd <- as.matrix(ipd_rig[,-1], nrow = 13)
p_ipd <-  matrix(1, nrow = 13, ncol = 1) %*%  as.matrix(ipd_rpg[1,-1], nrow = 1)
r_ipd <- r_ipd/p_ipd

r_gpd <- as.matrix(gpd_rig[,-1], nrow = 13)
p_gpd <- matrix(1, nrow = 13, ncol = 1) %*% as.matrix(gpd_rpg[1,-1], nrow = 1)
r_gpd <- r_gpd/p_gpd

diff_ipd_gpd <- r_ipd - r_gpd

diff_ipd_gpd2 <- t(diff_ipd_gpd[c(4,7,10,13),]) %>%
  as.data.frame()

names(diff_ipd_gpd2) <- c("h3", "h6", "h9", "h12")

### plot
# 1 quarter
p1<- ggplot(diff_ipd_gpd2, aes(x = h3)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_gpd2$h3 < 0), digits = 3))) +
  ggtitle("h = 3")+
  xlab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

# 2 quarter
p2 <- ggplot(diff_ipd_gpd2, aes(x = h6)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_gpd2$h6 < 0), digits = 3))) +
  ggtitle("h = 6")+
  xlab("")+
  ylab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

# 3 quarter
p3 <- ggplot(diff_ipd_gpd2, aes(x = h9)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_gpd2$h9 < 0), digits = 3))) +
  ggtitle("h = 9")+
  xlab("")+
  ylab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

# 4 quarter
p4 <- ggplot(diff_ipd_gpd2, aes(x = h12)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_gpd2$h12 < 0), digits = 3))) +
  ggtitle("h = 12")+
  xlab("")+
  ylab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

ggarrange(p1, p2, p3, p4, ncol = 4)
ggsave("modfig/irfs_ipd-gpd.png", device = "png", width = 9, height = 2.25)

##### difference reaction of rig to ipd and rpg shock
r_ipd <- as.matrix(ipd_rig[,-1], nrow = 13)
p_ipd <-  matrix(1, nrow = 13, ncol = 1) %*%  as.matrix(ipd_rpg[1,-1], nrow = 1)
r_ipd <- r_ipd/p_ipd

r_rpg <- as.matrix(rpg_rig[,-1], nrow = 13)
p_rpg <- matrix(1, nrow = 13, ncol = 1) %*% as.matrix(rpg_rpg[1,-1], nrow = 1)
r_rpg <- r_rpg/p_rpg

diff_ipd_rpg <- r_ipd - r_rpg

diff_ipd_rpg2 <- t(diff_ipd_rpg[c(4,7,10,13),]) %>%
  as.data.frame()

names(diff_ipd_rpg2) <- c("h3", "h6", "h9", "h12")

### plot
# 1 quarter
p1 <- ggplot(diff_ipd_rpg2, aes(x = h3)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_rpg2$h3 < 0), digits = 3))) +
  ggtitle("h = 3")+
  xlab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

# 2 quarter
p2 <- ggplot(diff_ipd_rpg2, aes(x = h6)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_rpg2$h6 < 0), digits = 3))) +
  ggtitle("h = 6")+
  xlab("")+
  ylab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

# 3 quarter
p3 <- ggplot(diff_ipd_rpg2, aes(x = h9)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_rpg2$h9 < 0), digits = 3))) +
  ggtitle("h = 9")+
  xlab("")+
  ylab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

# 4 quarter
p4 <- ggplot(diff_ipd_rpg2, aes(x = h12)) +
  geom_density(fill = "dark grey", col = "black", alpha = 0.6) +
  geom_vline(aes(xintercept = 0), col = "red")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(caption = paste("p = ",round(mean(diff_ipd_rpg2$h12 < 0), digits = 3))) +
  ggtitle("h = 12")+
  xlab("")+
  ylab("")+
  xlim(-2,2)+ 
  theme(text= element_text(family="LM Roman 10"))

ggarrange(p1, p2, p3, p4, ncol = 4)
ggsave("modfig/irfs_ipd-rpg.png", device = "png", width = 9, height = 2.25)

##################################################################################################
###################################### Robustness Checks #########################################
##################################################################################################

nc <- 4    # number of robustness checks
al <- 0.32 # confidence level for bounds

###### generate data frame for results

### price elasticity of drilling
elast_rob <- data.frame(model = c("Baseline", "Lower a34", "Higher a34", "Higher elasticity", "No trend"),
                        p11_low = rep(0,nc+1), p11_med = rep(0,nc+1), p11_high = rep(0,nc+1), p11_p = rep(0,nc+1),
                        p12_low = rep(0,nc+1), p12_med = rep(0,nc+1), p12_high = rep(0,nc+1), p12_p = rep(0,nc+1),
                        p13_low = rep(0,nc+1), p13_med = rep(0,nc+1), p13_high = rep(0,nc+1), p13_p = rep(0,nc+1),
                        p14_low = rep(0,nc+1), p14_med = rep(0,nc+1), p14_high = rep(0,nc+1), p14_p = rep(0,nc+1),
                        p21_low = rep(0,nc+1), p21_med = rep(0,nc+1), p21_high = rep(0,nc+1), p21_p = rep(0,nc+1),
                        p22_low = rep(0,nc+1), p22_med = rep(0,nc+1), p22_high = rep(0,nc+1), p22_p = rep(0,nc+1),
                        p23_low = rep(0,nc+1), p23_med = rep(0,nc+1), p23_high = rep(0,nc+1), p23_p = rep(0,nc+1),
                        p24_low = rep(0,nc+1), p24_med = rep(0,nc+1), p24_high = rep(0,nc+1), p24_p = rep(0,nc+1))

# add data from baseline model
elast_rob[1,2:4]   <- quantile(diff_ipd_gpd2$h3, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,5]     <- round(mean(diff_ipd_gpd2$h3 < 0), digits = 3)
elast_rob[1,6:8]   <- quantile(diff_ipd_gpd2$h6, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,9]     <- round(mean(diff_ipd_gpd2$h6 < 0), digits = 3)
elast_rob[1,10:12] <- quantile(diff_ipd_gpd2$h9, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,13]    <- round(mean(diff_ipd_gpd2$h9 < 0), digits = 3)
elast_rob[1,14:16] <- quantile(diff_ipd_gpd2$h12, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,17]    <- round(mean(diff_ipd_gpd2$h12 < 0), digits = 3)
elast_rob[1,18:20] <- quantile(diff_ipd_rpg2$h3, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,21]    <- round(mean(diff_ipd_rpg2$h3 < 0), digits = 3)
elast_rob[1,22:24] <- quantile(diff_ipd_rpg2$h6, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,25]    <- round(mean(diff_ipd_rpg2$h6 < 0), digits = 3)
elast_rob[1,26:28] <- quantile(diff_ipd_rpg2$h9, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,29]    <- round(mean(diff_ipd_rpg2$h9 < 0), digits = 3)
elast_rob[1,30:32] <- quantile(diff_ipd_rpg2$h12, probs = c(al/2, 0.5, 1-al/2))
elast_rob[1,33]    <- round(mean(diff_ipd_rpg2$h12 < 0), digits = 3)

### demand side contribution
contrib_rob <- data.frame(model = c("Baseline", "Lower a34", "Higher a34", "Higher elasticity", "No trend"),
                          p1_low = rep(0,nc+1), p1_med = rep(0,nc+1), p1_high = rep(0,nc+1),
                          p2_low = rep(0,nc+1), p2_med = rep(0,nc+1), p2_high = rep(0,nc+1),
                          p3_low = rep(0,nc+1), p3_med = rep(0,nc+1), p3_high = rep(0,nc+1))

# add data from baseline model
contrib_rob[1,2:4]   <- quantile(win1, probs = c(al/2, 0.5, 1-al/2))
contrib_rob[1,5:7]   <- quantile(win2, probs = c(al/2, 0.5, 1-al/2))
contrib_rob[1,8:10]  <- quantile(win3, probs = c(al/2, 0.5, 1-al/2))

###### bound for reaction of ipd to gas demand shock lower
rob_res <- bayes_4_4_var(star = as.Date("1993-11-01"), en = as.Date("2020-01-01"), M, K, a34 = 0.0027)
elast_rob[2,2:33] <- f_elast_rob(rob_res, al)
contrib_rob[2,2:10] <- f_contrib_rob(rob_res, 6, s1, e1, s2, e2, s3, e3, al)

###### bound for reaction of ipd to gas demand shock higher
rob_res <- bayes_4_4_var(star = as.Date("1993-11-01"), en = as.Date("2020-01-01"), M, K/3, a34 = 0.0053) # higher acceptance rate thus less rotations
elast_rob[3,2:33] <- f_elast_rob(rob_res, al)
contrib_rob[3,2:10] <- f_contrib_rob(rob_res, 6, s1, e1, s2, e2, s3, e3, al)

###### bound for price elasticity of natural gas production higher
rob_res <- bayes_4_4_var(star = as.Date("1993-11-01"), en = as.Date("2020-01-01"), M, K/2, ela = 0.085)
elast_rob[4,2:33] <- f_elast_rob(rob_res, al)
contrib_rob[4,2:10] <- f_contrib_rob(rob_res, 6, s1, e1, s2, e2, s3, e3, al)

###### no broken trend
rob_res <- bayes_4_4_var(star = as.Date("1993-11-01"), en = as.Date("2020-01-01"), M, K, tren = FALSE)
elast_rob[5,2:33] <- f_elast_rob(rob_res, al)
contrib_rob[5,2:10] <- f_contrib_rob(rob_res, 6, s1, e1, s2, e2, s3, e3, al)

############# plot results

### price elasticity of drilling

# Plot 1
elast_rob$model <- factor(elast_rob$model, levels = c("No trend", "Higher elasticity", "Higher a34", "Lower a34", "Baseline"))
p_pos <- 0.95
fonts <- 3.5
xl_mi <- -0.8
xl_ma <- 1.4

p11 <- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p11_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p11_low, xmax = p11_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p11_p),
         hjust = 0, size = fonts,
         fontface = ifelse(elast_rob$p11_p == elast_rob$p11_p[1], "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 3")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p12<- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p12_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p12_low, xmax = p12_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p12_p),
            hjust = 0, size = fonts,
            fontface = ifelse(elast_rob$p12_p == elast_rob$p12_p[1], "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 6")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p13 <- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p13_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p13_low, xmax = p13_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p13_p),
            hjust = 0, size = fonts,
            fontface = ifelse(elast_rob$p13_p == elast_rob$p13_p[1], "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 9")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p14 <- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p14_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p14_low, xmax = p14_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p14_p),
            hjust = 0, size = fonts,
            fontface = ifelse(elast_rob$p14_p == elast_rob$p14_p[1], "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 12")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

ggarrange(p11, p12, p13, p14, ncol = 4, widths = c(1.45, rep(1,4)))
ggsave("modfig/rob_irfs_ipd-gpd.png", device = "png", width = 9, height = 2.25)

# Plot 2
p_pos <- 1.3
fonts <- 3.5
xl_mi <- -0.15
xl_ma <- 1.7

p21 <- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p21_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p21_low, xmax = p21_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p21_p),
            hjust = 0, size = fonts,
            fontface = ifelse(seq(1,5) == 1, "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 3")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p22<- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p22_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p22_low, xmax = p22_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p22_p),
            hjust = 0, size = fonts,
            fontface = ifelse(seq(1,5) == 1, "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 6")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p23 <- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p23_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p23_low, xmax = p23_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p23_p),
            hjust = 0, size = fonts,
            fontface = ifelse(seq(1,5) == 1, "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 9")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p24 <- ggplot(elast_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p24_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p24_low, xmax = p24_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 0), col = "red")+
  geom_text(aes(x = p_pos, label = p24_p),
            hjust = 0, size = fonts,
            fontface = ifelse(seq(1,5) == 1, "bold", "plain")) +
  xlab("") +
  ylab("") +
  ggtitle("h = 12")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

ggarrange(p21, p22, p23, p24, ncol = 4, widths = c(1.45, rep(1,4)))
ggsave("modfig/rob_irfs_ipd-rpg.png", device = "png", width = 9, height = 2.25)

### demand side contributions 
contrib_rob$model <- factor(elast_rob$model, levels = c("No trend", "Higher elasticity", "Higher a34", "Lower a34", "Baseline"))
xl_mi <- -20
xl_ma <- 110

p1 <- ggplot(contrib_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p1_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p1_low, xmax = p1_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 50), col = "red", linewidth = 0.5) +
  geom_vline(aes(xintercept = 75), col = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("Percent") +
  ylab("") +
  ggtitle("2008M06-2009M09") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p2 <- ggplot(contrib_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p2_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p2_low, xmax = p2_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 50), col = "red", linewidth = 0.5) +
  geom_vline(aes(xintercept = 75), col = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("Percent") +
  ylab("") +
  ggtitle("2008M06-2012M04") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

p3 <- ggplot(contrib_rob, aes(y = model, size = model, linewidth = model)) +
  geom_point(aes(x = p3_med)) +
  scale_size_manual(values = c(rep(0.8,nc), 1.4)) +
  geom_linerange(aes(xmin = p3_low, xmax = p3_high)) +
  scale_linewidth_manual(values = c(rep(0.35,nc), 0.7)) +
  geom_vline(aes(xintercept = 50), col = "red", linewidth = 0.5) +
  geom_vline(aes(xintercept = 75), col = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("Percent") +
  ylab("") +
  ggtitle("2007M06-2013M06") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.y = element_blank()) +
  xlim(xl_mi, xl_ma) +
  theme(text= element_text(family="LM Roman 10"))

ggarrange(p1, p2, p3, ncol = 3, widths = c(1.43, 1,1))
ggsave("modfig/rob_contrib.png", device = "png", width = 7.5, height = 2.25)

save(elast_rob, file = "moddata/elast_rob.RData")
save(contrib_rob, file = "moddata/contrib_rob.RData")
