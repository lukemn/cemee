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
library(dae)
library(RColorBrewer)


# Results of the simulations can be loaded below
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

proj_M=eigen(gamma)$vectors

Ne=1000
df_trates=NULL
list_P_mat=list()

df_trates_drift=NULL
list_P_mat_drift=list()

df_trates_beta=NULL
list_P_mat_beta=list()

initial_P=matrix(rep(0,6),6,1)
G_to_simulate=VCV_mat[[3]]$G1_mat/2
G_to_simulate_rotated <- t(proj_M)%*% G_to_simulate%*%proj_M

sav_all_P=list()
sav_all_P_drift=list()
sav_all_P_beta=list()

sav_df_trates=list()
sav_df_trates_drift=list()
sav_df_trates_beta=list()

#beta_coef <- c(2,.75,.75,.75,.75,2)
#max_dist_opt <- c(0,.5,2,2,.5,.15)

rotated_phenotypes <- NULL
for(i in 1:nrow(final_fertility)) rotated_phenotypes <- rbind(rotated_phenotypes,t(t(proj_M)%*%t(final_fertility[i,9:14])))
dim(rotated_phenotypes)

max_dist_opt <- colSds(rotated_phenotypes)
beta_coef <-  3*(eigen(gamma)$values)
max_fit_diff <- 0.05

for(rep in 1:50){

current_P=initial_P
current_P_drift=initial_P
current_P_beta=initial_P

all_P=t(current_P)
all_P_drift=t(current_P)
all_P_beta=t(current_P)

k_list=1

exp_current_P=exp(current_P + sav_Means_P_values)
list_P_mat[[k_list]]=(t(matrix(c(-sum(exp_current_P[1:2]), exp_current_P[1:3],-sum(exp_current_P[3:4]), exp_current_P[4:6],-sum(exp_current_P[5:6])),3,3)))
df_trates=expm(20*list_P_mat[[k_list]])[1,]

exp_current_P_drift= exp_current_P
list_P_mat_drift[[k_list]]=(t(matrix(c(-sum(exp_current_P_drift[1:2]), exp_current_P_drift[1:3],-sum(exp_current_P_drift[3:4]), exp_current_P_drift[4:6],-sum(exp_current_P_drift[5:6])),3,3)))
df_trates_drift=expm(20*list_P_mat_drift[[k_list]])[1,]

exp_current_P_beta= exp_current_P
list_P_mat_beta[[k_list]]=(t(matrix(c(-sum(exp_current_P_beta[1:2]), exp_current_P_beta[1:3],-sum(exp_current_P_beta[3:4]), exp_current_P_beta[4:6],-sum(exp_current_P_beta[5:6])),3,3)))
df_trates_beta=expm(20*list_P_mat_beta[[k_list]])[1,]


for(i in 1:100){
k_list=k_list+1

# Drift and gamma
all_P=rbind(all_P, t(current_P)+rmvnorm((G_to_simulate)%*%(gamma)%*%((current_P)-(rep(0,6))),(G_to_simulate)/Ne))
current_P=matrix(all_P[nrow(all_P),],6,1)
exp_current_P=exp(current_P + sav_Means_P_values)

list_P_mat[[k_list]]=(t(matrix(c(-sum(exp_current_P[1:2]), exp_current_P[1:3],-sum(exp_current_P[3:4]), exp_current_P[4:6],-sum(exp_current_P[5:6])),3,3)))

if(i<=100) df_trates=rbind(df_trates,expm(10*list_P_mat[[k_list]])[1,])
if(i>100) df_trates=rbind(df_trates, df_trates[nrow(df_trates),])
# Drift only
all_P_drift=rbind(all_P_drift, t(current_P_drift)+rmvnorm(rep(0,6), G_to_simulate/Ne))
current_P_drift=matrix(all_P_drift[nrow(all_P_drift),],6,1)

exp_current_P_drift=exp(current_P_drift + sav_Means_P_values)
list_P_mat_drift[[k_list]]=(t(matrix(c(-sum(exp_current_P_drift[1:2]), exp_current_P_drift[1:3],-sum(exp_current_P_drift[3:4]), exp_current_P_drift[4:6],-sum(exp_current_P_drift[5:6])),3,3)))
df_trates_drift=rbind(df_trates_drift,expm(20*list_P_mat_drift[[k_list]])[1,])

# With beta and gamma

current_P_beta_rotated <- t(proj_M)%*%(matrix(current_P_beta,6,1))
is_beta_applied <- ifelse(abs(current_P_beta_rotated)>max_dist_opt,1,0)*sign(current_P_beta_rotated)
beta_vect =  is_beta_applied *(beta_coef * abs(current_P_beta_rotated)-max_dist_opt)

rot_select_matrix <- diag(eigen(gamma)$values) - beta_vect%*%t(beta_vect)
select_matrix <- proj_M%*% rot_select_matrix%*%t(proj_M)

all_P_beta=rbind(all_P_beta, t(current_P_beta)+rmvnorm(G_to_simulate%*% select_matrix%*% current_P_beta, G_to_simulate/Ne))
current_P_beta=matrix(all_P_beta[nrow(all_P_beta),],6,1)


exp_current_P_beta=exp(current_P_beta + sav_Means_P_values)
list_P_mat_beta[[k_list]]=(t(matrix(c(-sum(exp_current_P_beta[1:2]), exp_current_P_beta[1:3],-sum(exp_current_P_beta[3:4]), exp_current_P_beta[4:6],-sum(exp_current_P_beta[5:6])),3,3)))
df_trates_beta=rbind(df_trates_beta,expm(20*list_P_mat_beta[[k_list]])[1,])


}

sav_all_P[[rep]]=all_P
sav_all_P_drift[[rep]]=all_P_drift
sav_all_P_beta[[rep]]=all_P_beta

sav_df_trates[[rep]]=df_trates
sav_df_trates_drift[[rep]]=df_trates_drift
sav_df_trates_beta[[rep]]=df_trates_beta

}



for(i in 1:length(sav_all_P)){
	for(k in 1:nrow(sav_all_P[[i]])){
sav_all_P[[i]][k,]= sav_all_P[[i]][k,]+ sav_Means_P_values
sav_all_P_drift[[i]][k,]= sav_all_P_drift[[i]][k,]+ sav_Means_P_values
sav_all_P_beta[[i]][k,]= sav_all_P_beta[[i]][k,]+ sav_Means_P_values

	}
}

#INITIAL FREQS
exp_current_P_init=exp(sav_Means_P_values)
list_P_mat_init=(t(matrix(c(-sum(exp_current_P_init[1:2]), exp_current_P_init[1:3],-sum(exp_current_P_init[3:4]), exp_current_P_init[4:6],-sum(exp_current_P_init[5:6])),3,3)))
df_trates_init=expm(10* list_P_mat_init)[1,]

ylim_vect <- cbind(
c(-2.2,-3.3,-1.5,-5.5,-1.3,-3.4,-2.5,-1.5,-4),
c(-0.6,-1.2,0,-2.9,.8,-1,.5,2.5,-2))



##### Beta histogramms ?

vect_dist_opt <- NULL
vect_dist_opt_drift <- NULL
vect_dist_opt_beta <- NULL
for(k in 1:6){
temp_vect1 <-NULL;temp_vect2 <-NULL;temp_vect3 <-NULL
for(i in 1:length(sav_all_P)){
temp_vect1 <- c(temp_vect1,sav_all_P[[i]][101,k] - sav_Means_P_values[k])
temp_vect2 <- c(temp_vect2, sav_all_P_drift[[i]][101,k] - sav_Means_P_values[k])
temp_vect3 <- c(temp_vect3, sav_all_P_beta[[i]][101,k] - sav_Means_P_values[k])
}
vect_dist_opt <- rbind(vect_dist_opt, temp_vect1/100)
vect_dist_opt_drift <- rbind(vect_dist_opt_drift, temp_vect2/100)
vect_dist_opt_beta <- rbind(vect_dist_opt_beta, temp_vect3/100)

}



##################################################

sav_all_P_rotated= sav_all_P
sav_all_P_drift_rotated=sav_all_P_drift
sav_all_P_beta_rotated=sav_all_P_beta
for(i in 1:length(sav_all_P)){
	for(k in 1:nrow(sav_all_P[[i]])){
sav_all_P_rotated[[i]][k,]= t(proj_M)%*%matrix(sav_all_P_rotated[[i]][k,],6,1)
sav_all_P_drift_rotated[[i]][k,]= t(proj_M)%*%matrix(sav_all_P_drift_rotated[[i]][k,],6,1)
sav_all_P_beta_rotated[[i]][k,]= t(proj_M)%*%matrix(sav_all_P_beta_rotated[[i]][k,],6,1)

	}
}

sav_Means_P_values_rot = t(proj_M)%*% sav_Means_P_values

# Save here the results of the simulations for Figures
save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq_OU-Process.RData")


rm(list=ls())
load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq_OU-Process.RData')

#Fig 8A

par(mfrow=c(1,3),bty="n",yaxt="s",xaxt="s",mar=c(5,4,4,2))

k=1

plot(sav_df_trates[[1]][1:nrow(sav_df_trates[[1]]),k],type="n",ylim= c(.4,.8),ylab="Proportion of staying still",xlab="Generations",main='Drift + Gamma')
for(i in 1:50) lines(sav_df_trates[[i]][1:nrow(sav_df_trates[[1]]),k],type="l",col="black")

plot(sav_df_trates[[1]][1:nrow(sav_df_trates[[1]]),k],type="n",ylim= c(.4,.8),ylab="Proportion of staying still",xlab="Generations",main='Drift only')
for(i in 1:50) lines(sav_df_trates_drift[[i]][1:nrow(sav_df_trates_drift[[1]]),k],type="l",col="red")

plot(sav_df_trates[[1]][1:nrow(sav_df_trates[[1]]),k],type="n",ylim= c(.4,.8),ylab="Proportion of staying still",xlab="Generations",main='Drift + Gamma + Beta')
for(i in 1:50) lines(sav_df_trates_beta[[i]][1:nrow(sav_df_trates_beta[[1]]),k],type="l",col="green")


#Fig 8B
par(mfrow=c(1,3))
for(k in 1:3){
	if(k ==1){
	temp_dist_opt <- vect_dist_opt
	vect_main <- 'Drift + Gamma'	
	} 
	if(k ==2){
	 temp_dist_opt <- vect_dist_opt_drift
	 vect_main <- 'Drift only'	
	 }
	if(k ==3){
		temp_dist_opt <- vect_dist_opt_beta
		vect_main <- 'Drift + Gamma + Beta'	
		}
i=1
plot(density(temp_dist_opt[i,]),col=brewer.pal(n = 6, name = 'Dark2')[i],bty="n",las=1,xlab=expression(paste("Net directional selection gradient (", beta[net], ")")),xlim=c(-.01,0.01),ylab="Density",main= vect_main)
#mtext(side=3,"Transition rates")
for(i in 2:6){
	 lines(density(temp_dist_opt[i,]),col=brewer.pal(n = 6, name = 'Dark2')[i])
	 }
abline(v=rowMedians(temp_dist_opt),col=brewer.pal(n = 6, name = 'Dark2'))

if(k==1) legend(-.01,250,c("SF","SB","FS","FB","BS","BF"),lwd=2,col=brewer.pal(n = 6, name = 'Dark2'))
}

#Fig B9

par(mfrow=c(2,3))
new_ylim_vect=cbind(sav_Means_P_values_rot-.8,sav_Means_P_values_rot+.8)
vect_dist_opt = cbind(sav_Means_P_values_rot-max_dist_opt,sav_Means_P_values_rot+ max_dist_opt)


for(k in 1:6){

plot(sav_all_P_rotated[[1]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",ylim= new_ylim_vect[k,],ylab=expression(y[i]),xlab="Generations")
for(i in 2:50) lines(sav_all_P_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="black")
for(i in 1:50) lines(sav_all_P_drift_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="red")

#for(i in 1:50) lines(sav_all_P_beta_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="green")


#abline(h= sav_Means_P_values_rot[k],col="blue")
#abline(h= vect_dist_opt[k,],col="blue",lty=2)

#points(rep(0,nrow(rotated_phenotypes)),rotated_phenotypes[,k]+sav_Means_P_values_rot[k],pch=16,cex=.6)
mtext(paste0("Canonical axis ",1:6)[k],side=3)
}

## Figure Sup / TR
par(mfrow=c(2,3),bty="n",yaxt="s",xaxt="s",mar=c(5,4,4,2))
for(k in 1:6){

if(k %in% c(1,4)){
	 plot(sav_all_P[[1]][1:nrow(sav_all_P[[1]]),k],type="l",ylim= ylim_vect[k,],ylab="log transition rate",xlab="Generations")
	 }else{
	 plot(sav_all_P[[1]][1:nrow(sav_all_P[[1]]),k],type="l",ylim= ylim_vect[k,],ylab="",xlab="Generations")	 	
	 	}

for(i in 2:50) lines(sav_all_P[[i]][1:nrow(sav_all_P[[1]]),k],type="l",col="black")
for(i in 1:50) lines(sav_all_P_drift[[i]][1:nrow(sav_all_P[[1]]),k],type="l",col="red")
#for(i in 1:50) lines(sav_all_P_beta[[i]][1:nrow(sav_all_P[[1]]),k],type="l",col="green")

abline(h= sav_Means_P_values[k],col="blue")
mtext(c(
'SF',"SB","FS","FB","BS","BF")[k],side=3)
}
