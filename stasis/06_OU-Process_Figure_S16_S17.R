#quartz()
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

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")


Ne=1000
df_trates=NULL
list_P_mat=list()
df_trates_drift=NULL
list_P_mat_drift=list()
initial_P=matrix(rep(0,6),6,1)
G_to_simulate=VCV_mat[[1]]$G1_mat/2

sav_all_P=list()
sav_all_P_drift=list()
sav_df_trates=list()
sav_df_trates_drift=list()

for(rep in 1:50){

current_P=initial_P
current_P_drift=initial_P
all_P=t(current_P)
all_P_drift=t(current_P)

exp_current_P=exp(current_P + sav_Means_P_values)
#exp_current_P=exp(current_P)
k_list=1
list_P_mat[[k_list]]=(t(matrix(c(-sum(exp_current_P[1:2]), exp_current_P[1:3],-sum(exp_current_P[3:4]), exp_current_P[4:6],-sum(exp_current_P[5:6])),3,3)))

df_trates=expm(20*list_P_mat[[k_list]])[1,]

exp_current_P_drift= exp_current_P
list_P_mat_drift[[k_list]]=(t(matrix(c(-sum(exp_current_P_drift[1:2]), exp_current_P_drift[1:3],-sum(exp_current_P_drift[3:4]), exp_current_P_drift[4:6],-sum(exp_current_P_drift[5:6])),3,3)))
df_trates_drift=expm(20*list_P_mat_drift[[k_list]])[1,]


for(i in 1:100){
all_P=rbind(all_P, t(current_P)+rmvnorm((G_to_simulate)%*%(gamma)%*%((current_P)-(rep(0,6))),(G_to_simulate)/Ne))
current_P=matrix(all_P[nrow(all_P),],6,1)

exp_current_P=exp(current_P + sav_Means_P_values)

k_list=k_list+1
list_P_mat[[k_list]]=(t(matrix(c(-sum(exp_current_P[1:2]), exp_current_P[1:3],-sum(exp_current_P[3:4]), exp_current_P[4:6],-sum(exp_current_P[5:6])),3,3)))
df_trates=rbind(df_trates,expm(10*list_P_mat[[k_list]])[1,])

all_P_drift=rbind(all_P_drift, t(current_P_drift)+rmvnorm(rep(0,6), G_to_simulate/Ne))
current_P_drift=matrix(all_P_drift[nrow(all_P_drift),],6,1)

exp_current_P_drift=exp(current_P_drift + sav_Means_P_values)
list_P_mat_drift[[k_list]]=(t(matrix(c(-sum(exp_current_P_drift[1:2]), exp_current_P_drift[1:3],-sum(exp_current_P_drift[3:4]), exp_current_P_drift[4:6],-sum(exp_current_P_drift[5:6])),3,3)))
df_trates_drift=rbind(df_trates_drift,expm(20*list_P_mat_drift[[k_list]])[1,])

}

sav_all_P[[rep]]=all_P
sav_all_P_drift[[rep]]=all_P_drift
sav_df_trates[[rep]]=df_trates
sav_df_trates_drift[[rep]]=df_trates_drift

}

dim(df_trates_drift)

for(i in 1:length(sav_all_P)){
	for(k in 1:nrow(sav_all_P[[i]])){
sav_all_P[[i]][k,]= sav_all_P[[i]][k,]+ sav_Means_P_values
sav_all_P_drift[[i]][k,]= sav_all_P_drift[[i]][k,]+ sav_Means_P_values

	}
}

#INITIAL FREQS
exp_current_P_init=exp(sav_Means_P_values)
list_P_mat_init=(t(matrix(c(-sum(exp_current_P_init[1:2]), exp_current_P_init[1:3],-sum(exp_current_P_init[3:4]), exp_current_P_init[4:6],-sum(exp_current_P_init[5:6])),3,3)))
df_trates_init=expm(10* list_P_mat_init)[1,]

ylim_vect <- cbind(c(-2.2,-2.7,-2.4,-5.5,-1.3,-3.4,-2.5,-1.5,-4),
c(-0.6,-1.2,0,-2.9,.8,-1,.5,2.5,-2))

#Figures S17

par(mfrow=c(2,3),bty="n",yaxt="s",xaxt="s",mar=c(5,4,4,2))
for(k in 1:6){

if(k %in% c(1,4)){
	 plot(sav_all_P[[1]][1:nrow(sav_all_P[[1]]),k],type="l",ylim= ylim_vect[k,],ylab="log transition rate",xlab="Generations")
	 }else{
	 plot(sav_all_P[[1]][1:nrow(sav_all_P[[1]]),k],type="l",ylim= ylim_vect[k,],ylab="",xlab="Generations")	 	
	 	}

for(i in 2:50) lines(sav_all_P[[i]][1:nrow(sav_all_P[[1]]),k],type="l",col="black")
for(i in 1:50) lines(sav_all_P_drift[[i]][1:nrow(sav_all_P[[1]]),k],type="l",col="red")
abline(h= sav_Means_P_values[k],col="blue")
mtext(c(
'SF',"SB","FS","FB","BS","BF")[k],side=3)
}

sav_all_P_rotated= sav_all_P
sav_all_P_drift_rotated=sav_all_P_drift
for(i in 1:length(sav_all_P)){
	for(k in 1:nrow(sav_all_P[[i]])){
sav_all_P_rotated[[i]][k,]= t(proj_M)%*%matrix(sav_all_P_rotated[[i]][k,],6,1)
sav_all_P_drift_rotated[[i]][k,]= t(proj_M)%*%matrix(sav_all_P_drift_rotated[[i]][k,],6,1)
	}
}

sav_Means_P_values_rot = t(proj_M)%*% sav_Means_P_values


# Fig S16

ylim_vect_rot=t(proj_M)%*%ylim_vect[1:6,1]
ylim_vect_rot=cbind(ylim_vect_rot,t(proj_M)%*%ylim_vect[1:6,2])
ylim_vect_rot[2,]=c(-2,2)
par(mfrow=c(1,2),bty="n",yaxt="s",xaxt="s",mar=c(5,4,4,2))

k=1

plot(sav_all_P_rotated[[1]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",ylim= c(-.4,.1),ylab=expression(y[1]),xlab="Generations")
for(i in 2:50) lines(sav_all_P_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="black")
for(i in 1:50) lines(sav_all_P_drift_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="red")
abline(h= sav_Means_P_values_rot[k],col="blue")
mtext(paste0("Canonical axis ",1:6)[k],side=3)

k=6

plot(sav_all_P_drift_rotated[[1]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",ylim=c( .4,.7),ylab=expression(y[6]),yaxt="s",xlab="Generations")
#axis(side=2,at=c(.4,.6,.8,1))
for(i in 2:50) lines(sav_all_P_drift_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="red")
for(i in 1:50) lines(sav_all_P_rotated[[i]][1:nrow(sav_all_P_rotated[[1]]),k],type="l",col="black")
abline(h= sav_Means_P_values_rot[k],col="blue")
mtext(paste0("Canonical axis ",1:6)[k],side=3)
