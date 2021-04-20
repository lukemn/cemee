rm(list=ls())
gc()
#load libraries
library(rstan)
library(coda)
library(matrixStats)
library(data.table)
library(lme4)
library(MCMCglmm)
library(ggplot2)


load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData')
load("~/PATH/TO/DIR/Cemee_Pop_WI/WI_G-matrix-all.RData")

########## Projection of the G-matrices along the y1-6 #########



proj_M = eigen(gamma)$vectors


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



####################################

vect_col =c("black","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick","forestgreen")
vect_y_pops=rep(c(0,50,100,-33),c(1,3,3,1))

vect_margin=rbind(c(1.5,4,1.5,0),c(1.5,2,1.5,2),c(1.5,0,1.5,4),
c(2,4,1,0),c(2,2,1,2),c(2,0,1,4))
adjust_x=c(0,rep(c(-5,0,5),2),0)

vect_main=c(expression(y[1]),expression(y[2]),expression(y[3]),
expression(y[4]),expression(y[5]),expression(y[6]))
par(mfrow=c(2:3))
for(k in 1:6){

par(mar= vect_margin[k,])
temp_pop=subset(list_pts, Rotated_traits==k & !is.na(Variance))
post_modes_temp=tapply(temp_pop$Variance,temp_pop$pop,posterior.mode)


plot(post_modes_temp[c(1:7,10)]~vect_y_pops,ylim=c(0.002,0.3),type="n",log="y",xlim=c(-35,120),xlab="Generation",ylab="",bty="n",xaxt="n",yaxt="n",main= vect_main[k],cex.main=2)

if(k>3) axis(side=1,at=c(0,50,100),cex.axis=1.2)
if(k %in% c(1,4)){
 axis(side=2,at=c(0.005,0.01,0.05,0.1,0.3),las=1,cex.axis=1.2)
}else{
 axis(side=2,at=c(0.005,0.01,0.05,0.1,0.3),labels=rep("",5),las=1,cex.axis=1.2)
}


vect90=t(matrix(unlist(tapply(temp_pop$Variance,temp_pop$pop,function(x){
	HPDinterval(as.mcmc(x)) })),nrow=2))

vect80=t(matrix(unlist(tapply(temp_pop$Variance,temp_pop$pop,function(x){
	HPDinterval(as.mcmc(x),prob=.8) })),nrow=2))	
		
	arrows((vect_y_pops+adjust_x), vect80[,1], (vect_y_pops+adjust_x), vect80[,2],	code=0,length=0,angle=90,col= vect_col,lwd=3)

	arrows((vect_y_pops+adjust_x), vect90[,1], (vect_y_pops+adjust_x), vect90[,2],	code=3,length=.07,angle=90,col= grey(.6))

		
Ne_est=1000

ylimit_polygon= (1-2/(Ne_est))^c(0,100,100,0)*rep(vect80[1,],each=2)
polygon(x=c(0,100,100,0),y= ylimit_polygon,col=adjustcolor(grey(.2), alpha.f = 0.3),border=NA)



lines(0:100,(1-2/(Ne_est))^(0:100)*median(post_modes_temp[1]),lwd=2,lty=2)
points(post_modes_temp~I(vect_y_pops+ adjust_x),pch=16,col= vect_col)


}


save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/G_matrix_projection_for_plot.RData")

