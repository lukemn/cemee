#quartz()
rm(list = ls())
library(matrixStats)
gc()



OU_50_replicate <- function(empty){
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
load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData')

Ne = 1000
df_trates = NULL
list_P_mat = list()
df_trates_drift = NULL
list_P_mat_drift = list()

initial_P = matrix(rep(0, 6), 6, 1)

# A6140
G_to_simulate = VCV_mat[[1]]$G1_mat/2

sav_all_P = list()
sav_all_P_drift = list()
sav_df_trates = list()
sav_df_trates_drift = list()

for (rep in 1:50) {

	current_P = initial_P
	current_P_drift = initial_P
	all_P = t(current_P)
	all_P_drift = t(current_P)

	exp_current_P = exp(current_P + sav_Means_P_values)

	k_list = 1
	list_P_mat[[k_list]] = (t(matrix(c(-sum(exp_current_P[1:2]), exp_current_P[1:3], 
		-sum(exp_current_P[3:4]), exp_current_P[4:6], -sum(exp_current_P[5:6])), 
		3, 3)))

	df_trates = expm(20 * list_P_mat[[k_list]])[1, ]

	exp_current_P_drift = exp_current_P
	list_P_mat_drift[[k_list]] = (t(matrix(c(-sum(exp_current_P_drift[1:2]), 
		exp_current_P_drift[1:3], -sum(exp_current_P_drift[3:4]), exp_current_P_drift[4:6], 
		-sum(exp_current_P_drift[5:6])), 3, 3)))
	df_trates_drift = expm(20 * list_P_mat_drift[[k_list]])[1, ]


	for (i in 1:200) {
		all_P = rbind(all_P, t(current_P) + rmvnorm((G_to_simulate) %*% 
			(gamma) %*% ((current_P) - (rep(0, 6))), (G_to_simulate)/Ne))
		current_P = matrix(all_P[nrow(all_P), ], 6, 1)

		exp_current_P = exp(current_P + sav_Means_P_values)
		
		k_list = k_list + 1
		list_P_mat[[k_list]] = (t(matrix(c(-sum(exp_current_P[1:2]), 
			exp_current_P[1:3], -sum(exp_current_P[3:4]), exp_current_P[4:6], 
			-sum(exp_current_P[5:6])), 3, 3)))
		df_trates = rbind(df_trates, expm(10 * list_P_mat[[k_list]])[1, 
			])

		all_P_drift = rbind(all_P_drift, t(current_P_drift) + rmvnorm(rep(0, 
			6), G_to_simulate/Ne))
		current_P_drift = matrix(all_P_drift[nrow(all_P_drift), ], 6, 
			1)

		exp_current_P_drift = exp(current_P_drift + sav_Means_P_values)
		list_P_mat_drift[[k_list]] = (t(matrix(c(-sum(exp_current_P_drift[1:2]), 
			exp_current_P_drift[1:3], -sum(exp_current_P_drift[3:4]), 
			exp_current_P_drift[4:6], -sum(exp_current_P_drift[5:6])), 
			3, 3)))
		df_trates_drift = rbind(df_trates_drift, expm(20 * list_P_mat_drift[[k_list]])[1, 
			])

	}

	sav_all_P[[rep]] = all_P
	sav_all_P_drift[[rep]] = all_P_drift
	sav_df_trates[[rep]] = df_trates
	sav_df_trates_drift[[rep]] = df_trates_drift

}

#for (i in 1:length(sav_all_P)) {
#	for (k in 1:nrow(sav_all_P[[i]])) {
#		sav_all_P[[i]][k, ] = sav_all_P[[i]][k, ] + sav_Means_P_values
#		sav_all_P_drift[[i]][k, ] = sav_all_P_drift[[i]][k, ] + sav_Means_P_values
#	}
#}

return(list(sav_df_trates, sav_df_trates_drift))
}

library(parallel)
nbcores=3
clust <- makeCluster(nbcores)
param_list=list(); for(i in 1:10) param_list[[i]] <- c("")
output_list <- parLapply(clust, param_list , OU_50_replicate)
stopCluster(clust)

par(mar=c(5,6,4,2))
mean_A6=rep(0,201)
mean_A6_drift=rep(0,201)
for(ii in 1:10){
temp_trates= output_list[[ii]][[1]]
all_trates=NULL
for(i in 1:50) all_trates=rbind(all_trates, temp_trates[[i]][,1])
A6_var=colVars(all_trates)

temp_trates=output_list[[ii]][[2]]
all_trates=NULL
for(i in 1:50) all_trates=rbind(all_trates, temp_trates[[i]][,1])
A6_drift_var=colVars(all_trates)
if(ii==1){
plot(A6_var[1:100],type='l',bty="n",lwd=1,ylab="Phenotypic variance between \nreplicate populations",xlab="Generations",col="gray")	
}else{
lines(A6_var[1:100],lwd=1,col="gray")
	}
lines(A6_drift_var[1:100],col="lightgreen",lwd=1)

mean_A6= mean_A6+ A6_var
mean_A6_drift = mean_A6_drift + A6_drift_var
}

lines(mean_A6[1:100]/ii,lwd=2)
lines(mean_A6_drift[1:100]/ii,lwd=2,col="forestgreen")


save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/OU-process_01.RData")


