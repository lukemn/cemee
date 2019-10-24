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


## Produce Figure 8, S11 and S13

load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData')

MA_lines = new.env()
load('~/PATH/TO/DIR/MA_lines/M_matrices_estimates.RData', envir = MA_lines)
load("~/PATH/TO/DIR/Cemee_Pop_WI/WI_G-matrix-all.RData")

temp<-cbind(c(1:7),c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100"))
df_for_rot_BP <- NULL

for(i in 1:7){
A_temp= VCV_mat[[as.numeric(temp[i,1])]]$G1_mat/2
A_temp <- I(t(proj_M)%*%A_temp%*%proj_M)
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= temp[i,2],mean_VCV=colSums(abs(A_temp)),EV=1:6))
}

A_temp=MA_lines$VCV_mat[[1]]$G1_mat/5
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= "N2",mean_VCV=colSums(abs(A_temp)),EV=1:6))

A_temp=MA_lines$VCV_mat[[2]]$G1_mat/5
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= "PB",mean_VCV=colSums(abs(A_temp)),EV=1:6))

A_temp= WI_VCV_mat$G1_mat/2
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= "WI",mean_VCV=colSums(abs(A_temp)),EV=1:6))

A_temp= WI_VCV_mat_No_N2$G1_mat/2
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= "WI_noN2",mean_VCV=colSums(abs(A_temp)),EV=1:6))

A_temp= WI_VCV_mat_No_N2_with_CX$G1_mat/2
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_BP=rbind(df_for_rot_BP,data.frame(pop= "WI_noN2_CX",mean_VCV=colSums(abs(A_temp)),EV=1:6))


names(df_for_rot_BP)[2] = "Variance"
names(df_for_rot_BP)[3] = "Rotated_traits"

#distribution of slopes

# We have 10.000 MCMC samples, we will split them in 100 x 100
list_pts=NULL

# Compute all the EV
list_proj_M=list()
all_EV_df=NULL
for(i in 1:nrow(post_beta[[1]])){

temp_vect <-post_beta[[1]][i,2:22]
rdm_gamma=  matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)
list_proj_M[[i]]= eigen(rdm_gamma)$vectors
all_EV_df=rbind(all_EV_df,eigen(rdm_gamma)$values)
}

rdm_spl=matrix(sample(1:10000),100,100)
all_pts_df=NULL
	for(j in 1:7){
temp_pop=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")[j]
all_Variance=NULL

for(k in 1:100){
	temp_index=sample(1:length(list_proj_M),1)
	temp_proj_M= list_proj_M[[temp_index]]
temp_Variance=NULL

		for(r in 1:100){
	temp_G=t(temp_proj_M)%*%I(matrix(VCV_mat[[j]]$VCV_Mat[rdm_spl[k,r],1:36],6,6)/2)%*% temp_proj_M
temp_Variance=rbind(temp_Variance,colSums(abs(temp_G)))
}
all_pts_df = rbind(all_pts_df,data.frame(pop=temp_pop,Variance=as.numeric(apply(temp_Variance,2,posterior.mode)),Rotated_traits=1:6,niter=k, EV_gamma= all_EV_df[temp_index,]))
}
}

for(j in 1:5){
temp_pop=c("N2","PB306","WI","WI_noN2","WI_noN2_CX")[j]
all_Variance=NULL

for(k in 1:100){
	temp_index=sample(1:length(list_proj_M),1)
	temp_proj_M= list_proj_M[[temp_index]]
temp_Variance=NULL

		for(r in 1:100){
		if(temp_pop=="WI"){
			temp_G=t(temp_proj_M)%*%I(matrix(WI_VCV_mat$VCV_Mat[rdm_spl[k,r],1:36],6,6)/2)%*%temp_proj_M
			 }else if(temp_pop=="WI_noN2"){
			temp_G=t(temp_proj_M)%*%I(matrix(WI_VCV_mat_No_N2$VCV_Mat[rdm_spl[k,r],1:36],6,6)/2)%*%temp_proj_M			 	
			 }else if(temp_pop=="WI_noN2_CX"){
			temp_G=t(temp_proj_M)%*%I(matrix(WI_VCV_mat_No_N2_with_CX$VCV_Mat[rdm_spl[k,r],1:36],6,6)/2)%*%temp_proj_M			 	
			 }else{
			 temp_G=t(temp_proj_M)%*%I(matrix(MA_lines$VCV_mat[[j]]$VCV_Mat[rdm_spl[k,r],1:36],6,6)/5)%*% temp_proj_M
			 	
			 	}

temp_Variance=rbind(temp_Variance,colSums(abs(temp_G)))
}
all_pts_df =rbind(all_pts_df,data.frame(pop=temp_pop,Variance=as.numeric(apply(temp_Variance,2,posterior.mode)),Rotated_traits=1:6,niter=k, EV_gamma= all_EV_df[temp_index,]))
}
}

df_for_rot_BP$EV_gamma=eigen(gamma)$values

#####

list_pts= all_pts_df
df_for_rot_BP$pop=as.character(df_for_rot_BP$pop)
df_for_rot_BP$pop[df_for_rot_BP$pop=='PB']="PB306"
df_for_rot_BP$pop=factor(df_for_rot_BP$pop,levels=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100","N2","PB306","WI", "WI_noN2","WI_noN2_CX"))
df_for_rot_BP$pop2=c(1,2,2,2,3,3,3,4,5,6,7,8)[as.numeric(df_for_rot_BP$pop)]
list_pts=merge(list_pts,unique(df_for_rot_BP[,c("pop","pop2")]))
vcol=rep(c(grey(.5),"cornflowerblue", "firebrick3","darkgreen","magenta","yellow","cyan","black"),1)


save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_Canonical_Gamma_regression_withWI_noN2.RData")


###  Fig 8

temp2=subset(list_pts,substring(pop,1,2)%in%c("A6","CA"))
rdm_points =sample(1:nrow(temp2),nrow(temp2))

plot(abs(temp2$EV_gamma[rdm_points]),temp2$Variance[rdm_points],log="x",xlim=c(0.015,20),ylim=c(0.01, .3),pch=16,col=adjustcolor(rep(vcol,c(1,3,3,1,1,1,1,1))[as.numeric(temp2$pop)][rdm_points],alpha=.5),cex=.8,bty="n",ylab=c("Amount of genetic variance"),xlab=expression(paste("Strength of phenotypic selection ( |",lambda[i],"| )")))


for(k in 2:3){
	k_pop =c(1:3)
	 for(i in k_pop){
i_pop=i
temp_pop=subset(list_pts,pop2== i_pop & !is.na(Variance))
median_pop_x=tapply(temp_pop $EV_gamma, temp_pop $Rotated_traits, posterior.mode)
median_pop_y=tapply(temp_pop $Variance, temp_pop $Rotated_traits, posterior.mode)


if(k==1){
	}else if(k==3){
		points(abs(median_pop_x),median_pop_y,pch=21,cex=1.5,bg=vcol[i],col="black")
	}else{


	x_predict=exp(seq(log(0.05),log(15),length.out=30))

temp_mod=lm(median_pop_y~log(abs(median_pop_x)))
lines(c(0.02,20),predict(temp_mod,data.frame(median_pop_x=c(0.02,20))),lwd=3,col=vcol[i])


}}}

### Fig S11 and S13

vcol[6]="black"
vcol[7]="turquoise"
vcol[8]="yellow"

for(i_plot in 1:2){

par(mar=c(5,4,4,2))
if(i_plot==1) temp2=subset(list_pts,substring(pop,1,2)%in%c("WI","WI_noN2","WI_noN2_CX"))
if(i_plot==2) temp2=subset(list_pts,!substring(pop,1,2)%in%c("A6","CA","WI"))

rdm_points =sample(1:nrow(temp2),nrow(temp2))

if(i_plot==1){

plot(abs(temp2$EV_gamma[rdm_points]),temp2$Variance[rdm_points],log="x",xlim=c(0.015,20),ylim=c(0.01, .44),pch=16,col=adjustcolor(rep(vcol,c(1,3,3,1,1,1,1,1))[as.numeric(temp2$pop)][rdm_points],alpha=.5),cex=.8,bty="n",ylab=c("Amount of genetic variance"),xlab=expression(paste("Strength of phenotypic selection ( |",lambda[i],"| )")))

}else{
#par(mar=c(5,1,4,3))
plot(abs(temp2$EV_gamma[rdm_points]),temp2$Variance[rdm_points],log="x",xlim=c(0.015,20),ylim=c(0.01, .44),pch=16,col=adjustcolor(rep(vcol,c(1,3,3,1,1,2,1,1))[as.numeric(temp2$pop)][rdm_points],alpha=.5),cex=.8,bty="n",xlab=expression(paste("Strength of phenotypic selection ( |",lambda[i],"| )")),ylab=c("Amount of genetic variance"))
	
	}

for(k in 2:3){
	if(i_plot==1) k_pop =c(6,7,8)
	if(i_plot==2) k_pop =c(4:5)
	 for(i in k_pop){
i_pop=i#levels(as.factor(df_for_rot_BP$pop2))[i]
temp_pop=subset(list_pts,pop2== i_pop & !is.na(Variance))
median_pop_x=tapply(temp_pop $EV_gamma, temp_pop $Rotated_traits, posterior.mode)
median_pop_y=tapply(temp_pop $Variance, temp_pop $Rotated_traits, posterior.mode)


if(k==1){

	}else if(k==3){
		points(abs(median_pop_x),median_pop_y,pch=21,cex=2,bg=vcol[i],col="white")
	}else{
	temp_pop=subset(list_pts,pop2== i_pop & !is.na(Variance))
	x_predict=exp(seq(log(0.05),log(15),length.out=30))

temp_pop$X=((abs(temp_pop$EV_gamma)))
temp_pop$Y=((abs(temp_pop$Variance)))


temp_mod=lm(median_pop_y~log(abs(median_pop_x)))
lines(c(0.02,20),predict(temp_mod,data.frame(median_pop_x=c(0.02,20))),lwd=3.1,col="white")
lines(c(0.02,20),predict(temp_mod,data.frame(median_pop_x=c(0.02,20))),lwd=3,col=vcol[i])


}
}

}


}





