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
library(ks)
library(gridExtra)
library(cowplot)


load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData')
proj_M=eigen(gamma)$vectors

MA_lines = new.env()
load('~/PATH/TO/DIR/MA_lines/M_matrices_estimates.RData', envir = MA_lines)
load("~/PATH/TO/DIR/Cemee_Pop_WI/WI_G-matrix-all.RData")

temp<-cbind(c(1:7),c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100"))
df_for_rot_BP <- NULL

for(i in 1:7){
A_temp= VCV_mat[[as.numeric(temp[i,1])]]$G1_mat/2
A_temp <- I(t(proj_M)%*%A_temp%*%proj_M)
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= temp[i,2],mean_VCV=diag(A_temp),EV=1:6))
}

A_temp= WI_VCV_mat$G1_mat/2
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= "WI",mean_VCV=diag(A_temp),EV=1:6))

names(df_for_rot_BP)[2] = "Variance"
names(df_for_rot_BP)[3] = "Rotated_traits"

#distribution of slopes

post_all=model_MCMC$Sol[,2:22]

# We have 1.500 MCMC samples, we will split them in 1 x 1500
list_pts=NULL

# Compute all the EV
list_proj_M=list()
all_EV_df=NULL
for(i in 1:nrow(post_all)){

temp_vect <-post_all[i,]
temp_vect= temp_vect[c(7:21,1:6)]
rdm_gamma=  matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

list_proj_M[[i]]= eigen(rdm_gamma)$vectors
all_EV_df=rbind(all_EV_df,eigen(rdm_gamma)$values)
}

rdm_spl=NULL
for(i in 1:400) rdm_spl=rbind(rdm_spl ,sample(1:1500,10))

all_pts_df=NULL
	for(j in 1:7){
temp_pop=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")[j]
all_Variance=NULL

for(k in 1:nrow(rdm_spl)){
	temp_index=sample(1:length(list_proj_M),1)
	temp_proj_M= list_proj_M[[temp_index]]
temp_Variance=NULL

		for(r in 1:ncol(rdm_spl)){
	temp_G=t(temp_proj_M)%*%I(matrix(VCV_mat[[j]]$VCV_Mat[rdm_spl[k,r],1:36],6,6)/2)%*% temp_proj_M
temp_Variance=rbind(temp_Variance,diag(temp_G))
}
all_pts_df = rbind(all_pts_df,data.frame(pop=temp_pop,Variance=as.numeric(apply(temp_Variance,2,posterior.mode)),Rotated_traits=1:6,niter=k, EV_gamma= all_EV_df[temp_index,]))
}
}


temp_pop="WI"
all_Variance=NULL

for(k in 1:nrow(rdm_spl)){
	temp_index=sample(1:length(list_proj_M),1)
	temp_proj_M= list_proj_M[[temp_index]]
temp_Variance=NULL

		for(r in 1:ncol(rdm_spl)){

			temp_G=t(temp_proj_M)%*%I(matrix(WI_VCV_mat$VCV_Mat[rdm_spl[k,r],1:36],6,6)/2)%*%temp_proj_M
			temp_Variance=rbind(temp_Variance,diag(temp_G))
}
all_pts_df =rbind(all_pts_df,data.frame(pop=temp_pop,Variance=as.numeric(apply(temp_Variance,2,posterior.mode)),Rotated_traits=1:6,niter=k, EV_gamma= all_EV_df[temp_index,]))
}

df_for_rot_BP$EV_gamma=eigen(gamma)$values

#####

list_pts= all_pts_df
df_for_rot_BP$pop=as.character(df_for_rot_BP$pop)
df_for_rot_BP$pop[df_for_rot_BP$pop=='PB']="PB306"
df_for_rot_BP$pop=factor(df_for_rot_BP$pop,levels=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100","N2","PB306","WI"))
df_for_rot_BP$pop2=c(1,2,2,2,3,3,3,4,5,6,7,8)[as.numeric(df_for_rot_BP$pop)]
list_pts=merge(list_pts,unique(df_for_rot_BP[,c("pop","pop2")]))
vcol=rep(c(grey(.5),"cornflowerblue", "firebrick3","darkgreen","magenta","yellow","cyan","black"),1)

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_Canonical_Gamma_regression.RData")

