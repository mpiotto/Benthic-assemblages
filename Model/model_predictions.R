


library(tidyr)
library(doBy)
library(tidyverse)
library(raster)
library(sfsmisc)
library(ggplot2)
library(rstan)
library(ggpubr)
library(mvtnorm)
library(matrixcalc)
library(SimDesign)
library(DHARMa)
library("rstan")
library(matrixStats)
library(tibble)
library(viridis)

rm(list = ls())

calc_f = function(post)
{
  f = rep(NA, ncol(post))
  for(i in 1:ncol(post))
  {
    tmp1 = density(post[,i], from = min(post[,i])*1.1, to = max(post[,i]*1.1))
    if(min(post[,i])*1.1 > 0)
    {tmp1 = density(post[,i], from = 0, to = max(post[,i]*1.1))}
    dd = data.frame(x = tmp1$x, y = tmp1$y)
    if(mean(post[,i]) > 0)
    {
      f[i] = integrate.xy(dd$x, dd$y, a = 0, b = max(dd$x))
    }else{
      if(max(dd$x) < 0)
      {
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = max(dd$x))
      }else{
        f[i] = integrate.xy(dd$x, dd$y, a = min(dd$x), b = 0)
      }
    }
  }
  return(f)
}

load("inputdata.RData")
load("model_fit.RData")

fit_summary <- summary(fit, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
draws <- rstan::extract(fit)


# Genus specific responses ----------------------------------------------

thetas = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", # occupancy  
           "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8") # cover|present

# b1: intercept (bathy: ~ 15, 1994, presence); B1 intercept (bathy: ~ 15, 1994, cover)
# b2: 98 (presence), B2: 98 (cover)  diferencial de intercepto
# b3: 09 (presence), B3: 09 (cover)  diferencial de intercepto
# b4: 20 (presence), B4: 20 (cover)  diferencial de intercepto
# b5: slope-bathy (presence), B5: slope-bathy (cover)
# b6: s98 (presence), B6: s98 (cover)  diferencial de pendiente
# b7: s09 (presence), B7: s09 (cover)  diferencial de pendiente
# b8: s20 (presence), B8: s20 (cover)  diferencial de pendiente


bs <- fit_summary[grepl("betas", rownames(fit_summary)),] 

ids = seq(from = 1, to = nrow(bs), by = ngen)
length(ids) == length(thetas)

names_ = c("94_psi", "98_psi","09_psi", "20_psi", "deep_psi", "d98_psi","d09_psi", "d20_psi",
           "94__mu", "98_mu", "09_mu", "20_mu", "deep_mu", "d98_mu", "d09_mu", "d20_mu")

#Per beta, obtain posteriors
post_betas = vector("list", length(thetas))
B = vector("list",length(thetas))
for(i in 1:length(thetas))
{
  init = ids[i]
  fin = ids[i]+(ngen-1)
  B[[i]] = bs[init:fin,]
  post_betas[[i]] = draws$betas[, init:fin]
}
names(B) = names_
names(post_betas) = names_

## Species-specific responses to depth (occupancy)
d1994 = post_betas[["deep_psi"]]
d1998 = post_betas[["deep_psi"]] + post_betas[["d98_psi"]]
d2009 = post_betas[["deep_psi"]] + post_betas[["d09_psi"]]
d2021 = post_betas[["deep_psi"]] + post_betas[["d20_psi"]]

mean94 = apply(d1994, 2, mean) 
L94 = apply(d1994, 2, function(x) quantile(x, probs = 0.025))
U94 = apply(d1994, 2, function(x) quantile(x, probs = 0.975))

mean98 = apply(d1998, 2, mean) 
L98 = apply(d1998, 2, function(x) quantile(x, probs = 0.025))
U98 = apply(d1998, 2, function(x) quantile(x, probs = 0.975))

mean09 = apply(d2009, 2, mean) 
L09 = apply(d2009, 2, function(x) quantile(x, probs = 0.025))
U09 = apply(d2009, 2, function(x) quantile(x, probs = 0.975))

mean21 = apply(d2021, 2, mean) 
L21 = apply(d2021, 2, function(x) quantile(x, probs = 0.025))
U21 = apply(d2021, 2, function(x) quantile(x, probs = 0.975))

b_deep <- data.frame(genus = genus,
                     group = group,
                     mean = c(mean94, mean98, mean09, mean21),
                     L = c(L94, L98, L09, L21),
                     U = c(U94, U98, U09, U21),
                     year = rep(c("1994", "1998", "2009", "2021"), each = ngen))

tmp = subset(b_deep, group %in% c("G2", "G3"))
tmp$group <- factor(tmp$group, levels = c("G2", "G3"))

p1=ggplot(tmp, aes(x = genus, y = mean, color = as.factor(year)))+
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.6), 
                 alpha = 0.65,
                 size = 3) +
  geom_point(size = 4.5, 
             alpha = 0.7, 
             position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "darkgrey", size=1)+
  ylab("Species responses to depth") +
  xlab("Species") +
  scale_colour_manual(values = c(viridis(4)))+
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black", size = 18, angle = 90), 
        axis.text.y = element_text(colour = "black", size = 18),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.position = "bottom",
        panel.spacing.y = unit(0.5, "cm"),
        panel.grid.minor = element_blank())

tmp = subset(b_deep, group %in% c("G3"))

Genus_B = vector("list",ngen)

for(i in 1:ngen)
{
  out = cbind(post_betas[[1]][,i], post_betas[[2]][,i], post_betas[[3]][,i], 
              post_betas[[4]][,i], post_betas[[5]][,i], post_betas[[6]][,i], 
              post_betas[[7]][,i], post_betas[[8]][,i],post_betas[[9]][,i],
              post_betas[[10]][,i], post_betas[[11]][,i], post_betas[[12]][,i],
              post_betas[[13]][,i], post_betas[[14]][,i], post_betas[[15]][,i],
              post_betas[[16]][,i])
  colnames(out) = thetas
  Genus_B[[i]] = out
}
names(Genus_B) = genus

niter = nrow(Genus_B[[1]])


############## Predictions #################--

Preds = vector("list", ngen)

draws_z <- qnorm(p = ppoints(500))
sigma_ranef = sqrt(draws$sigma_tr^2 + draws$sigma_tr_sp^2)

for(i in 1:ngen){
  
  tmp_b = Genus_B[[i]]
  bathy = seq(min(X[,5]), max(X[,5]), length.out = 10)
  #n_tran = 100
  Psi94 = matrix(NA, 10, niter)
  Psi98 = matrix(NA, 10, niter)
  Psi09 = matrix(NA, 10, niter)
  Psi21 = matrix(NA, 10, niter)
  mu94 = matrix(NA, 10, niter)
  mu98 = matrix(NA, 10, niter)
  mu09 = matrix(NA, 10, niter)
  mu21 = matrix(NA, 10, niter)
  
  cat(sprintf("Processing genus %d of %d...\n", i, ngen))
  
  for(j in 1:niter){
    
    draw_scaled = draws_z * sigma_ranef[j] 
    
    for(b in 1:length(bathy)){
      
      #cat(sprintf("Processing iter %d of %d...\n", j, niter))
      
       Psi94[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,5] * bathy[b] + draw_scaled))
       Psi98[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,2] + (tmp_b[j,5] + tmp_b[j,6]) * bathy[b] + draw_scaled))
       Psi09[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,3] + (tmp_b[j,5] + tmp_b[j,7]) * bathy[b] + draw_scaled))
       Psi21[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,4] + (tmp_b[j,5] + tmp_b[j,8]) * bathy[b] + draw_scaled))
    
    }
    
    mu94[,j] = plogis(tmp_b[j,9] + tmp_b[j,13] * bathy)
    mu98[,j] = plogis(tmp_b[j,9] + tmp_b[j,10] + (tmp_b[j,13] + tmp_b[j,14]) * bathy)
    mu09[,j] = plogis(tmp_b[j,9] + tmp_b[j,11] + (tmp_b[j,13] + tmp_b[j,15]) * bathy)
    mu21[,j] = plogis(tmp_b[j,9] + tmp_b[j,12] + (tmp_b[j,13] + tmp_b[j,16]) * bathy)
    
  }
  
  cov94 = Psi94 * mu94
  cov98 = Psi98 * mu98
  cov09 = Psi09 * mu09
  cov21 = Psi21 * mu21 

  preds = list(Psi94, Psi98, Psi09, Psi21, 
               mu94, mu98, mu09, mu21, 
               cov94, cov98, cov09, cov21)
  
  Preds[[i]] = preds
  
}

names(Preds)= genus

save(Preds, file = "Preds.RData")  
n = nrow(Preds[[1]][[1]])


################################# OCCUPANCY ############################--

out1 = vector("list", ngen)

for (i in 1:ngen) {
  occ94 = Preds[[i]][[1]]
  occ98 = Preds[[i]][[2]]
  occ09 = Preds[[i]][[3]]
  occ21 = Preds[[i]][[4]]
  matrix = matrix(NA, 4*n, 9)
  matrix[,1] = rep(genus[i], n * 4)
  matrix[,2] = rep(morpho[morph_id[i]], n * 4)
  matrix[,3] = rep(type[[i]], n * 4)
  matrix[,4] = rep("occupancy", n * 4)
  matrix[,5] = rep(c("1994", "1998", "2009", "2021"), each = n)
  matrix[,6] = c(rowMeans(occ94), rowMeans(occ98), rowMeans(occ09), rowMeans(occ21))
  matrix[,7] = c(rowQuantiles(occ94, probs = c(0.025)), rowQuantiles(occ98, probs = c(0.025)),
                 rowQuantiles(occ09, probs = c(0.025)), rowQuantiles(occ21, probs = c(0.025)))
  matrix[,8] = c(rowQuantiles(occ94, probs = c(0.975)), rowQuantiles(occ98, probs = c(0.975)),
                 rowQuantiles(occ09, probs = c(0.975)), rowQuantiles(occ21, probs = c(0.975)))
  matrix[,9] = rep(seq(min(X[,5]), max(X[,5]), length.out = n), 4)
  
  out1[[i]] = matrix
}

OUT1 <- as.data.frame(do.call(rbind, out1), stringsAsFactors = FALSE)
names(OUT1) <- c("genus", "morpho", "type","param", "year", "mean", "L", "U", "bathy")


########################### COVER|PRESENT ############################--

out2 = vector("list", ngen)

for (i in 1:ngen) {
  cov94 = Preds[[i]][[5]]
  cov98 = Preds[[i]][[6]]
  cov09 = Preds[[i]][[7]]
  cov21 = Preds[[i]][[8]]
  matrix = matrix(NA, n * 4, 9)
  matrix[,1] = rep(genus[i], n * 4)
  matrix[,2] = rep(morpho[morph_id[i]], n * 4)
  matrix[,3] = rep(type[[i]], n * 4)
  matrix[,4] = rep("coverp", n * 4)
  matrix[,5] = rep(c("1994", "1998", "2009", "2021"), each = n)
  matrix[,6] = c(rowMeans(cov94), rowMeans(cov98), rowMeans(cov09), rowMeans(cov21))
  matrix[,7] = c(rowQuantiles(cov94, probs = c(0.025)), rowQuantiles(cov98, probs = c(0.025)),
                 rowQuantiles(cov09, probs = c(0.025)), rowQuantiles(cov21, probs = c(0.025)))
  matrix[,8] = c(rowQuantiles(cov94, probs = c(0.975)), rowQuantiles(cov98, probs = c(0.975)),
                 rowQuantiles(cov09, probs = c(0.975)), rowQuantiles(cov21, probs = c(0.975)))
  matrix[,9] = rep(seq(min(X[,5]), max(X[,5]), length.out = n), 4)
  
  out2[[i]] = matrix
}

OUT2 <- as.data.frame(do.call(rbind, out2), stringsAsFactors = FALSE)
names(OUT2) <- c("genus", "morpho", "type","param", "year", "mean", "L", "U", "bathy")


out1$sbathy <- as.numeric(out1$bathy) * sd(Env$bathy) + mean(Env$bathy)
out2$sbathy <- as.numeric(out2$bathy) * sd(Env$bathy) + mean(Env$bathy)

com = subset(out1, type == "common")

com <- com %>%
  mutate(genus = recode(genus, "MOLGULA" = "Molgula spp.", "MALACOBELEMNON" = "Malacobelemnon",
                        "LATERNULA"="Laternula sp.", "SEROLIS" ="Paraserolis sp."))

ord = c("Molgula spp.", "Malacobelemnon", "Laternula sp.", "Paraserolis sp.")
com = com %>%
  mutate(genus = factor(genus, levels = ord))

Y_pres <- Y
Y_pres[(Y > 0)] <- 1
tmp_pres <- pivot_longer(as.data.frame(Y_pres), cols = 1:ncol(Y_pres), 
                         names_to = "genus", values_to = "pres")

Env$year[Env$year == "2020"] <- "2022"
Env$year[Env$year == "2009"] <- "2010"


obs_pres <- data.frame(genus = tmp_pres$genus,
                       pres = tmp_pres$pres,
                       depth = Env$depth,
                       year = as.numeric(Env$year)) %>%
  group_by(genus, depth, year) %>%
  summarise(freq = mean(pres, na.rm = TRUE), .groups = "drop") %>%
  mutate(genus = recode(genus,
                        "MOLGULA" = "Molgula spp.",
                        "MALACOBELEMNON" = "Malacobelemnon",
                        "LATERNULA" = "Laternula sp.",
                        "SEROLIS" = "Paraserolis sp."),
         genus = factor(genus, levels = ord))

obs_pres <- obs_pres %>%
  mutate(genus = recode(genus, "MOLGULA" = "Molgula spp.", "MALACOBELEMNON" = "Malacobelemnon",
                        "LATERNULA"="Laternula sp.", "SEROLIS" ="Paraserolis sp."))

obs_pres$depth <- gsub("m", "", obs_pres$depth)

p1=ggplot(com, aes(x = sbathy, y = mean, col = as.factor(year)))+
  theme_bw()+
  scale_colour_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  geom_ribbon(aes(ymin = L, ymax = U, fill = as.factor(year)), alpha = 0.3)+
  geom_line(size = 1)+
  ylab("Occypancy")+
  xlab("Depth (m)")+
 # geom_point(data = obs_pres %>%
 #              dplyr::semi_join(com, by = c("genus", "year")),
 #            aes(x = as.numeric(depth), y = freq), 
#             shape = 1, color = "black", size = 3, alpha = 0.7, inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  facet_grid(rows=vars(genus), cols=vars(year))+
  theme(axis.text.x = element_text(colour = "black", size = 18), 
        axis.text.y = element_text(colour = "black", size = 18),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.position = "bottom",
        panel.spacing.y = unit(0.5, "cm"),
        panel.grid.minor = element_blank())


asc = subset(out1, genus %in% c("CNEMIDOCARPA", "ASCIDIA", "CORELLA"))
asc <- asc %>%
  mutate(genus = recode(genus, "CNEMIDOCARPA" = "Cnemidocarpa", "ASCIDIA" = "Ascidia sp.",
                        "CORELLA" ="Corella sp."))

p1=ggplot(asc, aes(x = sbathy, y = mean, col = as.factor(year)))+
  theme_bw()+
  scale_colour_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  geom_ribbon(aes(ymin = L, ymax = U, fill = as.factor(year)), alpha = 0.3)+
  geom_line(size = 1)+
  ylab("Occypancy")+
  xlab("Depth (m)")+
  # geom_point(data = obs_pres %>%
  #              dplyr::semi_join(com, by = c("genus", "year")),
  #            aes(x = as.numeric(depth), y = freq), 
  #             shape = 1, color = "black", size = 3, alpha = 0.7, inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  facet_grid(rows=vars(genus), cols=vars(year))+
  theme(axis.text.x = element_text(colour = "black", size = 18), 
        axis.text.y = element_text(colour = "black", size = 18),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.position = "bottom",
        panel.spacing.y = unit(0.5, "cm"),
        panel.grid.minor = element_blank())


# Discrete depth estimates ---------------------------------------------------

thetas = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", # occupancy  
           "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8") # cover|present


bs <- fit_summary[grepl("betas", rownames(fit_summary)),] 

ids = seq(from = 1, to = nrow(bs), by = ngen)
length(ids) == length(thetas)#Length of covariates 

names_ = c("94_psi", "98_psi","09_psi", "20_psi", "deep_psi", "d98_psi","d09_psi", "d20_psi",
           "94__mu", "98_mu", "09_mu", "20_mu", "deep_mu", "d98_mu", "d09_mu", "d20_mu")

#Per beta obtain posteriors
post_betas = vector("list", length(thetas))
B = vector("list",length(thetas))
for(i in 1:length(thetas))
{
  init = ids[i]
  fin = ids[i]+(ngen-1)
  B[[i]] = bs[init:fin,] # una lista por beta, spp en cols
  post_betas[[i]] = draws$betas[, init:fin]
}
names(B) = names_
names(post_betas) = names_

#Make a post_beta per genus c(b1, b2, b3, b4, b5, b6, b7, b8, B1,B2,B3,B4,B5,B6,B7,B8)

Genus_B = vector("list",ngen)

for(i in 1:ngen)
{
  out = cbind(post_betas[[1]][,i], post_betas[[2]][,i], post_betas[[3]][,i], 
              post_betas[[4]][,i], post_betas[[5]][,i], post_betas[[6]][,i], 
              post_betas[[7]][,i], post_betas[[8]][,i],post_betas[[9]][,i],
              post_betas[[10]][,i], post_betas[[11]][,i], post_betas[[12]][,i],
              post_betas[[13]][,i], post_betas[[14]][,i], post_betas[[15]][,i],
              post_betas[[16]][,i])
  colnames(out) = thetas
  Genus_B[[i]] = out
}
names(Genus_B) = genus

niter = nrow(Genus_B[[1]])
Preds = vector("list", ngen)

nT = max(transect_id)

draws_z <- qnorm(p = ppoints(500))

sigma_ranef = sqrt(draws$sigma_tr^2 + draws$sigma_tr_sp^2)


Preds_hist = vector("list", ngen)

for(i in 1:ngen)
{
  tmp_b = Genus_B[[i]]
  bathy = (c(15, 20, 25, 30) - mean(Env$bathy))/sd(Env$bathy)
  Psi94 = matrix(NA, 4, niter)
  Psi98 = matrix(NA, 4, niter)
  Psi09 = matrix(NA, 4, niter)
  Psi21 = matrix(NA, 4, niter)
  mu94 = matrix(NA, 4, niter)
  mu98 = matrix(NA, 4, niter)
  mu09 = matrix(NA, 4, niter)
  mu21 = matrix(NA, 4, niter)
  
  cat(sprintf("Processing genus %d of %d...\n", i, ngen))
  
  for(j in 1:niter){
    
    draw_scaled = draws_z * sigma_ranef[j]
    
    for(b in 1:length(bathy)){
      
      Psi94[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,5] * bathy[b] + draw_scaled))
      Psi98[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,2] + (tmp_b[j,5] + tmp_b[j,6]) * bathy[b] + draw_scaled))
      Psi09[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,3] + (tmp_b[j,5] + tmp_b[j,7]) * bathy[b] + draw_scaled))
      Psi21[b,j] = mean(plogis(tmp_b[j,1] + tmp_b[j,4] + (tmp_b[j,5] + tmp_b[j,8]) * bathy[b] + draw_scaled))
      
    }
    
    mu94[,j] = plogis(tmp_b[j,9] + tmp_b[j,13] * bathy) 
    mu98[,j] = plogis(tmp_b[j,9] + tmp_b[j,10] + (tmp_b[j,13] + tmp_b[j,14]) * bathy)
    mu09[,j] = plogis(tmp_b[j,9] + tmp_b[j,11] + (tmp_b[j,13] + tmp_b[j,15]) * bathy)
    mu21[,j] = plogis(tmp_b[j,9] + tmp_b[j,12] + (tmp_b[j,13] + tmp_b[j,16]) * bathy)
    
  }
  
  cov94 = Psi94 * mu94
  cov98 = Psi98 * mu98
  cov09 = Psi09 * mu09
  cov21 = Psi21 * mu21 
  
  preds = list(Psi94, Psi98, Psi09, Psi21, 
               mu94, mu98, mu09, mu21, 
               cov94, cov98, cov09, cov21)
  
  Preds_hist[[i]] = preds
  
}

names(Preds_hist)= genus

n = nrow(Preds_hist[[1]][[1]])

save(Preds_hist,
     file = "Preds_hist.RData")



## Dominant species occupancies ----------------------------------------------

##### MOLGULA

## 1994
mean(Preds_hist[["MOLGULA"]][[1]][1,])## occurrence 1994 at 15
mean(Preds_hist[["MOLGULA"]][[1]][2,])## occurrence 1994 at 20

mean(Preds_hist[["MOLGULA"]][[1]][3,]) ## occurrence 1994 at 25
quantile(Preds_hist[["MOLGULA"]][[1]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[1]][4,]) ## occurrence 1994 at 30
quantile(Preds_hist[["MOLGULA"]][[1]][4,], probs = c(0.025, 0.975))

tmp = colMeans(rbind(Preds_hist[["MOLGULA"]][[5]][3,], Preds_hist[["MOLGULA"]][[5]][4,])) ## cover|present 1994 at 25-30
mean(tmp)
quantile(tmp, probs = c(0.025, 0.975))

## 1998
mean(Preds_hist[["MOLGULA"]][[2]][3,]) ## occurrence 1998 at 25
quantile(Preds_hist[["MOLGULA"]][[2]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[2]][4,]) ## occurrence 1998 at 30
quantile(Preds_hist[["MOLGULA"]][[2]][4,], probs = c(0.025, 0.975))

tmp = rbind(Preds_hist[["MOLGULA"]][[1]][3,]/Preds_hist[["MOLGULA"]][[2]][3,],
            Preds_hist[["MOLGULA"]][[1]][4,]/Preds_hist[["MOLGULA"]][[2]][4,])
mean(colMeans(tmp)) ## comparing occupancy 1994/1998, 25-30 averaging
quantile(colMeans(tmp), probs = c(0.025, 0.975))


mean(Preds_hist[["MOLGULA"]][[6]][3,]) ## cover|p 1998 at 25
quantile(Preds_hist[["MOLGULA"]][[6]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[6]][4,]) ## cover|p 1998 at 30
quantile(Preds_hist[["MOLGULA"]][[6]][4,], probs = c(0.025, 0.975))

tmp = rbind(Preds_hist[["MOLGULA"]][[5]][3,]/Preds_hist[["MOLGULA"]][[6]][3,],
            Preds_hist[["MOLGULA"]][[5]][4,]/Preds_hist[["MOLGULA"]][[6]][4,])
mean(colMeans(tmp)) ## comparing cover|p 1994/1998 averaging 25-30m (CKKKKK)
quantile(colMeans(tmp), probs = c(0.025, 0.975))

## 2009
mean(Preds_hist[["MOLGULA"]][[3]][3,]) ## occupancy 2009 at 25
quantile(Preds_hist[["MOLGULA"]][[3]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[3]][4,])## occupancy 2009 at 30
quantile(Preds_hist[["MOLGULA"]][[3]][4,], probs = c(0.025, 0.975))

tmp = rbind(Preds_hist[["MOLGULA"]][[1]][3,]/Preds_hist[["MOLGULA"]][[3]][3,],
            Preds_hist[["MOLGULA"]][[1]][4,]/Preds_hist[["MOLGULA"]][[3]][4,])
mean(colMeans(tmp)) ## comparing occupancy 1994/2009
quantile(colMeans(tmp), probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[7]][3,]) ## cover|present 2009 at 25
quantile(Preds_hist[["MOLGULA"]][[7]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[7]][4,])## cover|present 2009 at 30
quantile(Preds_hist[["MOLGULA"]][[7]][4,], probs = c(0.025, 0.975))

tmp = rbind(Preds_hist[["MOLGULA"]][[5]][3,]/Preds_hist[["MOLGULA"]][[7]][3,],
            Preds_hist[["MOLGULA"]][[5]][4,]/Preds_hist[["MOLGULA"]][[7]][4,])
mean(colMeans(tmp)) ## comparing cover|present 1994/2009
quantile(colMeans(tmp), probs = c(0.025, 0.975))

## 2021
mean(Preds_hist[["MOLGULA"]][[4]][3,]) ## occupancy 2021 at 25
quantile(Preds_hist[["MOLGULA"]][[4]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[4]][4,])## occupancy 2021 at 30
quantile(Preds_hist[["MOLGULA"]][[4]][4,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[8]][3,]) ## cover|present 2021 at 25
quantile(Preds_hist[["MOLGULA"]][[8]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MOLGULA"]][[8]][4,])## cover|present 2021 at 30
quantile(Preds_hist[["MOLGULA"]][[8]][4,], probs = c(0.025, 0.975))


#### MALACOBELEMNON

## 1994
mean(Preds_hist[["MALACOBELEMNON"]][[1]][1,])## occurrence 1994 at 15
quantile(Preds_hist[["MALACOBELEMNON"]][[1]][1,], probs = c(0.025, 0.975))

mean(Preds_hist[["MALACOBELEMNON"]][[1]][2,])## occurrence 1994 at 20
quantile(Preds_hist[["MALACOBELEMNON"]][[1]][2,], probs = c(0.025, 0.975))

mean(Preds_hist[["MALACOBELEMNON"]][[1]][3,]) ## occurrence 1994 at 25
quantile(Preds_hist[["MALACOBELEMNON"]][[1]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MALACOBELEMNON"]][[1]][4,]) ## occurrence 1994 at 30
quantile(Preds_hist[["MALACOBELEMNON"]][[1]][4,], probs = c(0.025, 0.975))


## 1998
mean(Preds_hist[["MALACOBELEMNON"]][[2]][1,]) ## occurrence 1998 at 15
quantile(Preds_hist[["MALACOBELEMNON"]][[2]][1,], probs = c(0.025, 0.975))
mean(Preds_hist[["MALACOBELEMNON"]][[2]][2,]) ## occurrence 1998 at 20
mean(Preds_hist[["MALACOBELEMNON"]][[2]][3,]) ## occurrence 1998 at 25
mean(Preds_hist[["MALACOBELEMNON"]][[2]][4,]) ## occurrence 1998 at 30


## 2009
mean(Preds_hist[["MALACOBELEMNON"]][[3]][3,]) ## occupancy 2009 at 25
quantile(Preds_hist[["MALACOBELEMNON"]][[3]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MALACOBELEMNON"]][[3]][4,])## occupancy 2009 at 30
quantile(Preds_hist[["MALACOBELEMNON"]][[3]][4,], probs = c(0.025, 0.975))

## 2021
mean(Preds_hist[["MALACOBELEMNON"]][[4]][3,]) ## occupancy 2021 at 25
quantile(Preds_hist[["MALACOBELEMNON"]][[4]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["MALACOBELEMNON"]][[4]][4,])## occupancy 2021 at 30
quantile(Preds_hist[["MALACOBELEMNON"]][[4]][4,], probs = c(0.025, 0.975))


## average cover

tmp = rbind(Preds_hist[["MALACOBELEMNON"]][[5]][1,],
            Preds_hist[["MALACOBELEMNON"]][[5]][2,],
            Preds_hist[["MALACOBELEMNON"]][[5]][3,],
            Preds_hist[["MALACOBELEMNON"]][[5]][4,],
            Preds_hist[["MALACOBELEMNON"]][[6]][1,],
            Preds_hist[["MALACOBELEMNON"]][[6]][2,],
            Preds_hist[["MALACOBELEMNON"]][[6]][3,],
            Preds_hist[["MALACOBELEMNON"]][[6]][4,],
            Preds_hist[["MALACOBELEMNON"]][[7]][1,],
            Preds_hist[["MALACOBELEMNON"]][[7]][2,],
            Preds_hist[["MALACOBELEMNON"]][[7]][3,],
            Preds_hist[["MALACOBELEMNON"]][[7]][4,],
            Preds_hist[["MALACOBELEMNON"]][[8]][1,],
            Preds_hist[["MALACOBELEMNON"]][[8]][2,],
            Preds_hist[["MALACOBELEMNON"]][[8]][3,],
            Preds_hist[["MALACOBELEMNON"]][[8]][4,])

mean(colMeans(tmp))
quantile(tmp, probs = c(0.025, 0.975))

# Occupancy Molgula 1994 Malacobelemnon 2021 
tmp = colMeans(rbind(Preds_hist[["MALACOBELEMNON"]][[4]][3,]/Preds_hist[["MOLGULA"]][[4]][3,], #malaco 2021/molgula 1994 at 25
                     Preds_hist[["MALACOBELEMNON"]][[4]][4,]/Preds_hist[["MOLGULA"]][[4]][4,])) #malaco 2021/molgula 1994 at 30
mean(tmp)
quantile(tmp, probs = c(0.025, 0.975))


### LATERNULA

tmp=colMeans(rbind(Preds_hist[["LATERNULA"]][[2]][1,]/Preds_hist[["LATERNULA"]][[1]][1,],
                   Preds_hist[["LATERNULA"]][[2]][2,]/Preds_hist[["LATERNULA"]][[1]][2,]))## occurrence 1998/1994 at 15 and 20
mean(tmp)
quantile(tmp, probs = c(0.025, 0.975))

mean(Preds_hist[["LATERNULA"]][[4]][1,])## occurrence 2021 at 15
quantile(Preds_hist[["LATERNULA"]][[4]][1,],  probs = c(0.025, 0.975))

mean(Preds_hist[["LATERNULA"]][[4]][2,])## occurrence 2021 at 20
quantile(Preds_hist[["LATERNULA"]][[4]][2,], probs = c(0.025, 0.975))

mean(Preds_hist[["LATERNULA"]][[4]][3,])## occurrence 2021 at 25
quantile(Preds_hist[["LATERNULA"]][[4]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["LATERNULA"]][[4]][4,])## occurrence 2021 at 30
quantile(Preds_hist[["LATERNULA"]][[4]][4,], probs = c(0.025, 0.975))


## SEROLIS

tmp = rbind(Preds_hist[["SEROLIS"]][[3]][3,]/Preds_hist[["SEROLIS"]][[2]][3,],## 2009/1998 25 m
            Preds_hist[["SEROLIS"]][[3]][4,]/Preds_hist[["SEROLIS"]][[2]][4,])
mean(colMeans(tmp))
quantile(colMeans(tmp), probs = c(0.025, 0.975))

mean(Preds_hist[["SEROLIS"]][[4]][3,])## 2021 25 m
quantile(Preds_hist[["SEROLIS"]][[4]][3,], probs = c(0.025, 0.975))

mean(Preds_hist[["SEROLIS"]][[4]][4,])## 2021 30 m
quantile(Preds_hist[["SEROLIS"]][[4]][4,], probs = c(0.025, 0.975))


### CORELLA, ASCIDIA, CNEMIDOCARPA

tmp = rbind(Preds_hist[["CORELLA"]][[2]][4,]/Preds_hist[["CORELLA"]][[3]][4,],## 1998/2009 30 m
            Preds_hist[["ASCIDIA"]][[2]][4,]/Preds_hist[["ASCIDIA"]][[3]][4,])
            #Preds_hist[["CNEMIDOCARPA"]][[2]][4,]/Preds_hist[["CNEMIDOCARPA"]][[3]][4,])
mean(colMeans(tmp))
quantile(colMeans(tmp), probs = c(0.025, 0.975))

mean(Preds_hist[["CORELLA"]][[2]][4,]) ## Corella 1998
quantile(Preds_hist[["CORELLA"]][[2]][4,], probs = c(0.025, 0.975))
mean(Preds_hist[["CORELLA"]][[3]][4,]) # 2009
quantile(Preds_hist[["CORELLA"]][[3]][4,], probs = c(0.025, 0.975))
mean(Preds_hist[["CORELLA"]][[4]][4,]) # 2021
quantile(Preds_hist[["CORELLA"]][[4]][4,], probs = c(0.025, 0.975))

mean(Preds_hist[["CNEMIDOCARPA"]][[2]][4,]) ## Cnemidocarpa 1998
quantile(Preds_hist[["CNEMIDOCARPA"]][[2]][4,], probs = c(0.025, 0.975))
mean(Preds_hist[["CNEMIDOCARPA"]][[3]][4,]) # 2009
quantile(Preds_hist[["CNEMIDOCARPA"]][[3]][4,], probs = c(0.025, 0.975))
mean(Preds_hist[["CNEMIDOCARPA"]][[4]][4,]) # 2021
quantile(Preds_hist[["CNEMIDOCARPA"]][[4]][4,], probs = c(0.025, 0.975))

mean(Preds_hist[["ASCIDIA"]][[2]][4,]) ## Ascidia 1998
quantile(Preds_hist[["ASCIDIA"]][[2]][4,], probs = c(0.025, 0.975))
mean(Preds_hist[["ASCIDIA"]][[3]][4,]) # 2009
quantile(Preds_hist[["ASCIDIA"]][[3]][4,], probs = c(0.025, 0.975))
mean(Preds_hist[["ASCIDIA"]][[4]][4,]) # 2021
quantile(Preds_hist[["ASCIDIA"]][[4]][4,], probs = c(0.025, 0.975))

mean(Preds_hist[["CORELLA"]][[7]][4,]/Preds_hist[["CORELLA"]][[8]][4,]) # cover 2009/2021
quantile(Preds_hist[["CORELLA"]][[7]][4,]/Preds_hist[["CORELLA"]][[8]][4,], probs = c(0.025, 0.975))

mean(Preds_hist[["ASCIDIA"]][[7]][4,]/Preds_hist[["ASCIDIA"]][[8]][4,]) # cover 2009/2021
quantile(Preds_hist[["ASCIDIA"]][[7]][4,]/Preds_hist[["ASCIDIA"]][[8]][4,], probs = c(0.025, 0.975))

mean(Preds_hist[["CNEMIDOCARPA"]][[7]][4,]/Preds_hist[["CNEMIDOCARPA"]][[8]][4,]) # cover 2009/2021
quantile(Preds_hist[["CNEMIDOCARPA"]][[7]][4,]/Preds_hist[["CNEMIDOCARPA"]][[8]][4,], probs = c(0.025, 0.975))


# dominant ascidians mean occupancy in 1994 30m
ma94=colMeans(rbind(Preds_hist[["MOLGULA"]][[1]][4,], 
                    Preds_hist[["CORELLA"]][[1]][4,],
                    Preds_hist[["CNEMIDOCARPA"]][[1]][4,], 
                    Preds_hist[["ASCIDIA"]][[1]][4,]))

# dominant ascidians mean occupancy in 2021 30m
ma21=colMeans(rbind(Preds_hist[["MOLGULA"]][[4]][4,], 
                    Preds_hist[["CORELLA"]][[4]][4,],
                    Preds_hist[["CNEMIDOCARPA"]][[4]][4,], 
                    Preds_hist[["ASCIDIA"]][[4]][4,]))

mean(ma94/ma21)
quantile(ma94/ma21, probs = c(0.025, 0.975))

# dominant ascidians mean cove|p in 1994 30m
ma94=colMeans(rbind(Preds_hist[["MOLGULA"]][[5]][4,], 
                    Preds_hist[["CORELLA"]][[5]][4,],
                    Preds_hist[["CNEMIDOCARPA"]][[5]][4,], 
                    Preds_hist[["ASCIDIA"]][[5]][4,]))


# dominant ascidians mean cover|p in 2021 30m
ma21=colMeans(rbind(Preds_hist[["MOLGULA"]][[8]][4,], 
                    Preds_hist[["CORELLA"]][[8]][4,],
                    Preds_hist[["CNEMIDOCARPA"]][[8]][4,], 
                    Preds_hist[["ASCIDIA"]][[8]][4,]))

mean(ma94/ma21)
quantile(ma94/ma21, probs = c(0.025, 0.975))

rm(list = ls(pattern = "^tmp"))


## Cover plots

ls_coverp <- lapply(Preds_hist, function(sublista) {
  do.call(rbind, sublista[5:8])
})

coverp <- do.call(rbind, ls_coverp)

coverp <- data.frame(genus = rep(genus, each = nrow(ls_coverp[[1]])),
                     type = rep(type, each = nrow(ls_coverp[[1]])),
                     mean = rowMeans(coverp),
                     L = rowQuantiles(coverp, probs = 0.025),
                     U = rowQuantiles(coverp, probs = 0.975),
                     depth = rep(c(15, 20, 25, 30)),
                     year = rep(c("1994", "1998", "2009", "2021"), each = 4))

coverp = subset(coverp, !(year == "1998" & depth %in% c(15, 25)))
com = subset(coverp, type == "common" )

com <- com %>%
  mutate(genus = recode(genus, "MOLGULA" = "Molgula sp.", "MALACOBELEMNON" = "Malacobelemnon sp.",
                        "LATERNULA"="Laternula sp.", "SEROLIS" ="Paraserolis sp."))

Y_covp <- Y
Y_covp[Y == 0] <- NA
tmp_cp <- pivot_longer(as.data.frame(Y_covp), cols = 1:ncol(Y_covp), 
                         names_to = "genus", values_to = "covp")

Env$year[Env$year == "2020"] <- "2021"

obs_cp <- data.frame(genus = tmp_cp$genus,
                       covp = tmp_cp$covp,
                       depth = Env$depth,
                       year = as.numeric(Env$year))

obs_cp <- obs_cp %>%
  mutate(genus = recode(genus, "MOLGULA" = "Molgula sp.", "MALACOBELEMNON" = "Malacobelemnon sp.",
                        "LATERNULA"="Laternula sp.", "SEROLIS" ="Paraserolis sp."))

ord = c("Molgula sp.", "Malacobelemnon sp.", "Laternula sp.", "Paraserolis sp.")
com = com %>%
  mutate(genus = factor(genus, levels = ord))

obs_cp = obs_cp %>%
  mutate(genus = factor(genus, levels = ord))

obs_cp <- aggregate(covp ~ genus + depth + year, data = obs_cp, FUN = mean)
obs_cp$depth <- gsub("m", "", as.character(obs_cp$depth))

p1=ggplot(com, aes(x = as.factor(depth), y = mean, col = as.factor(year)))+
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.8,
                 size = 2) +
  geom_point(size = 4.5, 
             alpha = 0.9, 
             position = position_dodge(width = 0.5)) + 
 # geom_jitter(data = obs_cp, 
 #             aes(x = depth, y = covp, col = as.factor(year)), 
 #             shape = 1, size = 3, 
 #             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5))+
  scale_colour_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  ylab("Percentage conditional cover")+
  xlab("Depth (m)")+
  facet_grid(rows=vars(genus))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 13), 
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        strip.text = element_text(size=14),
        legend.position = "bottom",
        #panel.spacing.y = unit(0.5, "cm"),
        panel.grid.minor = element_blank())


asc = subset(coverp, genus %in% c("CNEMIDOCARPA", "ASCIDIA", "CORELLA"))
asc <- asc %>%
  mutate(genus = recode(genus, "CNEMIDOCARPA" = "Cnemidocarpa sp.", "ASCIDIA" = "Ascidia sp.",
                        "CORELLA" ="Corella sp."))


ord = c("Ascidia sp.", "Cnemidocarpa sp.", "Corella sp.")
asc = asc %>%
  mutate(genus = factor(genus, levels = ord))


p2=ggplot(asc, aes(x = as.factor(depth), y = mean, col = as.factor(year)))+
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.8,
                 size = 2) +
  geom_point(size = 4.5, 
             alpha = 0.9, 
             position = position_dodge(width = 0.5)) + 
  # geom_jitter(data = obs_cp, 
  #             aes(x = depth, y = covp, col = as.factor(year)), 
  #             shape = 1, size = 3, 
  #             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5))+
  scale_colour_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  ylab("Percentage conditional cover")+
  xlab("Depth (m)")+
  facet_grid(rows=vars(genus))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 13), 
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        strip.text = element_text(size=14),
        legend.position = "bottom",
        #panel.spacing.y = unit(0.5, "cm"),
        panel.grid.minor = element_blank())


# Zs (trait effects, gamma) -----------------------------------------------

target_tr = colnames(TT)
z <- fit_summary[grepl("z", rownames(fit_summary)),]
nt = length(target_tr)

tocc <- draws$z[, 1:40] 
tcov <- draws$z[, 41:80] 
niter = nrow(tocc)

tocc[, c(1,6,11,16,21,26,31,36)] 
idx = seq(1, ncol(tocc), by = 5)

int <- tocc[, idx]                
pred <- tocc[, idx + 1]           
oppor <- tocc[, idx + 2]          
rep <- tocc[, idx + 3]            
hab <- tocc[, idx + 4]            
colnames(rep) = c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")
colnames(hab) = c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")

fil <- int 
oppor <- int + oppor
pred <- int + pred
list_gamma_occ = list(fil, pred, oppor)

fed = tr$fed

names(list_gamma_occ) = unique(fed)
for (i in seq_along(list_gamma_occ)) {
  colnames(list_gamma_occ[[i]]) = c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")
}


### Cover
int_cov <- tcov[, idx]                
pred_cov <- tcov[, idx + 1]           
oppor_cov <- tcov[, idx + 2]          
rep_cov <- tcov[, idx + 3]
hab_cov <- tcov[, idx + 4]
colnames(rep_cov) = c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")
colnames(hab_cov) = c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")
fil_cov <- int_cov # efecto de ser filtrador, sobre 94, dy98, dy09, dy21, bathy, s98, s09, s21
oppor_cov <- int_cov + oppor_cov
pred_cov <- int_cov + pred_cov
list_gamma_cov = list(fil_cov, pred_cov, oppor_cov)

fed = tr$fed

names(list_gamma_cov) = unique(fed)
for (i in seq_along(list_gamma_cov)) {
  colnames(list_gamma_cov[[i]]) = c("y94", "dy98", "dy09", "dy21", "d94", "s98", "s09", "s21")
}



# Habitat trait -----------------------------------------------------------

thetas = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", # occupancy  
           "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8") # cover|present

bs <- fit_summary[grepl("betas", rownames(fit_summary)),] 
ids = seq(from = 1, to = nrow(bs), by = ngen)
length(ids) == length(thetas)#Length of covariates  

names_ = c("94_psi", "98_psi","09_psi", "20_psi", "deep_psi", "d98_psi","d09_psi", "d20_psi",
           "94__mu", "98_mu", "09_mu", "20_mu", "deep_mu", "d98_mu", "d09_mu", "d20_mu")


#Per beta obtain posteriors
post_betas = vector("list", length(thetas))
B = vector("list",length(thetas))
for(i in 1:length(thetas))
{
  init = ids[i]
  fin = ids[i]+(ngen-1)
  B[[i]] = bs[init:fin,]
  post_betas[[i]] = draws$betas[, init:fin]
}
names(B) = names_
names(post_betas) = names_

hab_occ94 <- post_betas[["deep_psi"]]
hab_occ98 <- post_betas[["deep_psi"]]+post_betas[["d98_psi"]]
hab_occ09 <- post_betas[["deep_psi"]]+post_betas[["d09_psi"]]
hab_occ20 <- post_betas[["deep_psi"]]+post_betas[["d20_psi"]]

hab_occ <- data.frame(mean = c(colMeans(hab_occ94), colMeans(hab_occ98), 
                               colMeans(hab_occ09), colMeans(hab_occ20)),
                      L = c(colQuantiles(hab_occ94, probs = c(0.025)), colQuantiles(hab_occ98, probs = c(0.025)),
                            colQuantiles(hab_occ09, probs = c(0.025)), colQuantiles(hab_occ20, probs = c(0.025))),
                      U = c(colQuantiles(hab_occ94, probs = c(0.975)), colQuantiles(hab_occ98, probs = c(0.975)),
                            colQuantiles(hab_occ09, probs = c(0.975)), colQuantiles(hab_occ20, probs = c(0.975))),
                      hab = rep(TT[,5] * sd(TT[,5]) + mean(TT[,5]), length.out = length(genus)*4),
                      year = rep(c("1994", "1998", "2009", "2022"), each = length(genus)),
                      genus = rep(genus, length.out = length(genus)*4),
                      morfo = rep(tr$morfo, length.out = length(genus)*4))

hab_occ$shab <- hab_occ$hab * sd(tr$habitat) + mean(tr$habitat)

p=ggplot(hab_occ, aes(x = shab,  y = mean, color = year)) +    
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.55,
                 size = 0.5) +
  geom_point(size = 4, 
             alpha = 0.7, 
             position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, size = 1.3, linetype = "dashed", color = "darkgrey")+
  scale_colour_manual(values = c(viridis(4)))+
  #scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  ylab("Response to depth (occupancy)")+
  xlab("Habitat use")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  #facet_grid(rows = vars(year))+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14),
        legend.position = "none")


hab_occ$shab <- hab_occ$hab * sd(tr$habitat) + mean(tr$habitat)


# Cover|presence

hab_cov94 <- post_betas[["deep_mu"]]
hab_cov98 <- post_betas[["deep_mu"]]+post_betas[["d98_mu"]]
hab_cov09 <- post_betas[["deep_mu"]]+post_betas[["d09_mu"]]
hab_cov20 <- post_betas[["deep_mu"]]+post_betas[["d20_mu"]]

hab_cov <- data.frame(mean = c(colMeans(hab_cov94), colMeans(hab_cov98), 
                               colMeans(hab_cov09), colMeans(hab_cov20)),
                      L = c(colQuantiles(hab_cov94, probs = c(0.025)), colQuantiles(hab_cov98, probs = c(0.025)),
                            colQuantiles(hab_cov09, probs = c(0.025)), colQuantiles(hab_cov20, probs = c(0.025))),
                      U = c(colQuantiles(hab_cov94, probs = c(0.975)), colQuantiles(hab_cov98, probs = c(0.975)),
                            colQuantiles(hab_cov09, probs = c(0.975)), colQuantiles(hab_cov20, probs = c(0.975))),
                      hab = rep(TT[,5] * sd(TT[,5]) + mean(TT[,5]), length.out = length(genus)*4),
                      year = rep(c("1994", "1998", "2009", "2022"), each = length(genus)),
                      genus = rep(genus, length.out = length(genus)*4),
                      morfo = rep(tr$morfo, length.out = length(genus)*4))

hab_cov$shab <- hab_cov$hab * sd(tr$habitat) + mean(tr$habitat)
p=ggplot(hab_cov, aes(x = shab,  y = mean, color = year)) +    
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.55,
                 size = 0.5) +
  geom_point(size = 4, 
             alpha = 0.7, 
             position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, size = 1.3, linetype = "dashed", color = "darkgrey")+
  scale_colour_manual(values = c(viridis(4)))+
  #scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  ylab("Response to depth (cover|present)")+
  xlab("Habitat use")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  #facet_grid(rows = vars(year))+
  theme(axis.text.x = element_text(colour = "black", size = 13), 
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14),
        legend.position = "none")

# Reproductive trait ------------------------------------------------------------

thetas = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", # occupancy  
           "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8") # cover|present



bs <- fit_summary[grepl("betas", rownames(fit_summary)),] 
ids = seq(from = 1, to = nrow(bs), by = ngen)
length(ids) == length(thetas)#Length of covariates  


names_ = c("94_psi", "98_psi","09_psi", "20_psi", "deep_psi", "d98_psi","d09_psi", "d20_psi",
          "94__mu", "98_mu", "09_mu", "20_mu", "deep_mu", "d98_mu", "d09_mu", "d20_mu")


#Per beta obtain posteriors
post_betas = vector("list", length(thetas))
B = vector("list",length(thetas))
for(i in 1:length(thetas))
{
  init = ids[i]
  fin = ids[i]+(ngen-1)
  B[[i]] = bs[init:fin,]
  post_betas[[i]] = draws$betas[, init:fin]
}
names(B) = names_
names(post_betas) = names_

rep_occ94 <- post_betas[["deep_psi"]]
rep_occ98 <- post_betas[["deep_psi"]]+post_betas[["d98_psi"]]
rep_occ09 <- post_betas[["deep_psi"]]+post_betas[["d09_psi"]]
rep_occ20 <- post_betas[["deep_psi"]]+post_betas[["d20_psi"]]

TT[18,4] <- 1.4919332

rep_occ <- data.frame(mean = c(colMeans(rep_occ94), colMeans(rep_occ98), 
                               colMeans(rep_occ09), colMeans(rep_occ20)),
                      L = c(colQuantiles(rep_occ94, probs = c(0.025)), colQuantiles(rep_occ98, probs = c(0.025)),
                            colQuantiles(rep_occ09, probs = c(0.025)), colQuantiles(rep_occ20, probs = c(0.025))),
                      U = c(colQuantiles(rep_occ94, probs = c(0.975)), colQuantiles(rep_occ98, probs = c(0.975)),
                            colQuantiles(rep_occ09, probs = c(0.975)), colQuantiles(rep_occ20, probs = c(0.975))),
                      rep = rep(TT[,4] * sd(TT[,4]) + mean(TT[,4]), length.out = length(genus)*4),
                      year = rep(c("1994", "1998", "2009", "2022"), each = length(genus)),
                      genus = rep(genus, length.out = length(genus)*4),
                      morfo = rep(tr$morfo, length.out = length(genus)*4))

rep_occ$srep <- rep_occ$rep * sd(tr$rstrategy) + mean(tr$rstrategy)

p=ggplot(rep_occ, aes(x = srep,  y = mean, color = year)) +    
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.55,
                 size = 0.5) +
  geom_point(size = 4, 
             alpha = 0.7, 
             position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, size = 1.3, linetype = "dashed", color = "darkgrey")+
  scale_colour_manual(values = c(viridis(4)))+
  #scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  ylab("Response to depth (occupancy)")+
  xlab("Reproductive strategy")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 14), 
        axis.text.y = element_text(colour = "black", size = 14),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14),
        legend.position = "none")

# Cover|presence

rep_covp94 <- post_betas[["deep_mu"]]
rep_covp98 <- post_betas[["deep_mu"]]+post_betas[["d98_mu"]]
rep_covp09 <- post_betas[["deep_mu"]]+post_betas[["d09_mu"]]
rep_covp20 <- post_betas[["deep_mu"]]+post_betas[["d20_mu"]]

rep_covp <- data.frame(mean = c(colMeans(rep_covp94), colMeans(rep_covp98), 
                               colMeans(rep_covp09), colMeans(rep_covp20)),
                      L = c(colQuantiles(rep_covp94, probs = c(0.025)), colQuantiles(rep_covp98, probs = c(0.025)),
                            colQuantiles(rep_covp09, probs = c(0.025)), colQuantiles(rep_covp20, probs = c(0.025))),
                      U = c(colQuantiles(rep_covp94, probs = c(0.975)), colQuantiles(rep_covp98, probs = c(0.975)),
                            colQuantiles(rep_covp09, probs = c(0.975)), colQuantiles(rep_covp20, probs = c(0.975))),
                      rep = rep(TT[,4] * sd(TT[,4]) + mean(TT[,4]), length.out = length(genus)*4),
                      year = rep(c("1994", "1998", "2009", "2022"), each = length(genus)),
                      genus = rep(genus, length.out = length(genus)*4),
                      morfo = rep(tr$morfo, length.out = length(genus)*4))

rep_covp$srep <- rep_covp$rep * sd(tr$rstrategy) + mean(tr$rstrategy)

p=ggplot(rep_covp, aes(x = srep,  y = mean, color = year)) +    
  geom_linerange(aes(ymin = L, ymax = U), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.55,
                 size = 0.5) +
  geom_point(size = 4, 
             alpha = 0.7, 
             position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, size = 1.3, linetype = "dashed", color = "darkgrey")+
  scale_colour_manual(values = c(viridis(4)))+
  #scale_fill_manual(values = c("#414487FF", "#22A384FF", "#7AD151FF", "#FDE725FF"))+
  ylab("Response to depth (cover|present)")+
  xlab("Reproductive strategy")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  #facet_grid(rows = vars(year))+
  theme(axis.text.x = element_text(colour = "black", size = 13), 
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14),
        legend.position = "none")

ggsave(filename = "rep-cov.pdf", 
       path = "C:\\Users\\Administrador\\Desktop",
       plot = p,   # el objeto
       device = "pdf",    # el formato
       width = 11, height = 8, units = "cm",  
       dpi = 600)


# Rho paramater -----------------------------------------------------------

hist(draws$rho, main = "Rho parameter", xlab = "Rho")
abline(v = mean(draws$rho), col = "#414487FF", lw = 5)
dev.copy(png,"Rho.png",width=10,height=8,units="in",res=300)
dev.off()

mean(draws$rho)
quantile(draws$rho, probs = c(0.025, 0.975))

tmp <- data.frame(rho = draws$rho)

p=ggplot(tmp, aes(x = rho)) +
  geom_density(color = "black", fill = "grey", alpha = 0.65, size = 0.4) +
  geom_vline(xintercept = mean(draws$rho), color = "#414487FF", size = 1.2, linetype = "dashed") +
  labs(x = "ρ", y = "Frequency") +
  theme_bw()+
  #facet_grid(rows = vars(year))+
  theme(axis.text.x = element_text(colour = "black", size = 12), 
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=14))

