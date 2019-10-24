rm(list=ls());gc()
library(MCMCglmm)
library(psych)
library(ggplot2)
library(dplyr)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(corrplot)
library(Rmisc)
library(nlme)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

## We need to estimate the G-matrix of the WI
WI_for_G_estimation <- subset(final_export,population2=="WI")

for(j in c('temperature',"rel_humidity","logD")){
WI_for_G_estimation[,j] <- (WI_for_G_estimation[,j]-mean(WI_for_G_estimation[,j]))/sd(WI_for_G_estimation[,j])
}

	phen.var = diag(nb_trait) * diag(var(subset(WI_for_G_estimation, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC_WI <- MCMCglmm(cbind(c(T12_pred_02, T13_pred_02, T21_pred_02, T23_pred_02, T31_pred_02, T32_pred_02)) ~ (logD +rel_humidity+ temperature)+ trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = WI_for_G_estimation, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000)

post_dist = posterior.mode(model_MCMC_WI$VCV)
WI_G_mat=matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)

	WI_VCV_mat=list(Population = "WI", G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC_WI$VCV)

###### We also need a new estimates of the G without N2
WI_for_G_estimation <- subset(final_export,population2=="WI")

### We also need to retrieve here the data for the CX12311 line
original_ds = read.table("~/PATH/TO/DIRCemee_Pop_WI/Final_merged_data_NGM.txt",h=TRUE,sep="\t")
 data_CX12311 = subset(original_ds,data_group_name=="B400" & pop_label=="CX12311L0")
# We need the initial correction factors
cor_factors_forWI=as.numeric(WI_for_G_estimation[1,12:17]-WI_for_G_estimation[1,24:29])
k=0
for(i in vect_P_traits){
k=k+1
data_CX12311=cbind(data_CX12311,data_CX12311[,i]-cor_factors_forWI[k])
names(data_CX12311)[ncol(data_CX12311)] = paste0(i,'_pred_02')
}

data_CX12311$population="CX12311"
data_CX12311$population2="WI"

WI_for_G_estimation=rbind(WI_for_G_estimation[,names(WI_for_G_estimation)%in%names(data_CX12311)], data_CX12311[names(data_CX12311)%in%names(WI_for_G_estimation)])

for(j in c('temperature',"rel_humidity","logD")){
WI_for_G_estimation[,j] <- (WI_for_G_estimation[,j]-mean(WI_for_G_estimation[,j]))/sd(WI_for_G_estimation[,j])
}

# Remove the WIs
WI_for_G_estimation_CX=subset(WI_for_G_estimation,!pop_label%in%c("N2ancL0","CB4507L0"))

# Remove the CX
WI_for_G_estimation=subset(WI_for_G_estimation,!pop_label%in%c("N2ancL0","CB4507L0","CX12311L0"))

## Without the CX12311
	phen.var = diag(nb_trait) * diag(var(subset(WI_for_G_estimation, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC_WI_No_N2 <- MCMCglmm(cbind(c(T12_pred_02, T13_pred_02, T21_pred_02, T23_pred_02, T31_pred_02, T32_pred_02)) ~ (logD +rel_humidity+ temperature)+ trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = WI_for_G_estimation, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000)

post_dist = posterior.mode(model_MCMC_WI_No_N2$VCV)
WI_G_mat_No_N2=matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)

	WI_VCV_mat_No_N2=list(Population = "WI", G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC_WI_No_N2 $VCV)

## With the CX12311
	phen.var = diag(nb_trait) * diag(var(subset(WI_for_G_estimation_CX, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC_WI_No_N2_with_CX <- MCMCglmm(cbind(c(T12_pred_02, T13_pred_02, T21_pred_02, T23_pred_02, T31_pred_02, T32_pred_02)) ~ (logD +rel_humidity+ temperature)+ trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = WI_for_G_estimation_CX, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000)

post_dist = posterior.mode(model_MCMC_WI_No_N2_with_CX$VCV)
WI_G_mat_No_N2_with_CX=matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)

	WI_VCV_mat_No_N2_with_CX=list(Population = "WI_noN2_CX", G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC_WI_No_N2_with_CX$VCV)

### End of the estimation of G for WIs

save(list=c("WI_VCV_mat","WI_VCV_mat_No_N2","WI_VCV_mat_No_N2_with_CX"),file="~/PATH/TO/DIR/Cemee_Pop_WI/WI_G-matrix-all.RData")


