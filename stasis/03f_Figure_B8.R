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

load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_Canonical_Gamma_regression.RData')

temp<-cbind(c(1:2),c("N2","PB306"))
df_for_rot_MA <- NULL


for(i in 1:2){
A_temp= MA_lines$VCV_mat[[i]]$G1_mat/2
A_temp <- t(proj_M)%*%A_temp%*%proj_M
df_for_rot_MA=rbind(df_for_rot_MA,data.frame(pop= temp[i,2],mean_VCV=diag(A_temp),EV=1:6))
}

names(df_for_rot_MA)[2] = "Variance"
names(df_for_rot_MA)[3] = "Rotated_traits"

#distribution of slopes


all_pts_df_MA=NULL
	for(j in 1:2){
	temp_pop=temp[j,2]
	all_Variance=NULL

		for(k in 1:nrow(rdm_spl)){
		temp_index=sample(1:length(list_proj_M),1)
		temp_proj_M= list_proj_M[[temp_index]]
		temp_Variance=NULL

		for(r in 1:ncol(rdm_spl)){

		temp_G=t(temp_proj_M)%*%I(matrix(MA_lines$VCV_mat[[j]]$VCV_Mat[rdm_spl[k,r],1:36],6,6)/5)%*% temp_proj_M
		temp_Variance=rbind(temp_Variance,diag(temp_G))
}
	all_pts_df_MA = rbind(all_pts_df_MA,data.frame(pop=temp_pop,Variance=as.numeric(apply(temp_Variance,2,posterior.mode)),Rotated_traits=1:6,niter=k, EV_gamma= all_EV_df[temp_index,]))
}
}


df_for_rot_MA$EV_gamma=eigen(gamma)$values


list_pts_MA= all_pts_df_MA
df_for_rot_MA$pop=as.character(df_for_rot_MA$pop)

df_for_rot_MA$pop2=c(7,8)[as.numeric(as.factor(df_for_rot_MA$pop))]
list_pts_MA=merge(list_pts_MA,unique(df_for_rot_MA[,c("pop","pop2")]))



vect_col =c("black","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick","forestgreen")
vect_y_pops=rep(c(0,50,100,-33),c(1,3,3,1))

vect_margin=rbind(c(1.5,3,1.5,0),c(1.5,1.5,1.5,1.5),c(1.5,0,1.5,3),
c(2,3,1,0),c(2,1.5,1,1.5),c(2,0,1,3))
adjust_x=c(0,rep(c(-5,0,5),2),0)


vect_main=c(expression(y[1]),expression(y[2]),expression(y[3]),
expression(y[4]),expression(y[5]),expression(y[6]))
par(mfrow=c(1,1))
par(mar= c(5,4,4,2))

plot(post_modes_temp[c(1:7,10)]~vect_y_pops,ylim=c(0.,0.2),type="n",log="",xlim=c(0,7),xlab="Canonical Axis",ylab="Genetic variance",bty="n",xaxt="n",)
axis(side=1,at=1:6,labels=expression(y[1], y[2],y[3],y[4],y[5],y[6]))

for(k in 1:6){


	temp_pop=subset(list_pts, Rotated_traits==k & !is.na(Variance) & pop=="WI")
	post_modes_temp=tapply(temp_pop$Variance,temp_pop$pop,posterior.mode)

	temp_MA=subset(list_pts_MA, Rotated_traits==k & !is.na(Variance))
	post_modes_temp_MA=tapply(temp_MA$Variance,temp_MA$pop,posterior.mode)

vect90=t(matrix(unlist(tapply(temp_pop$Variance,temp_pop$pop,function(x){
	sort(x)[c(0.05,0.95)*length(x)]})),nrow=2))

vect90_MA=t(matrix(unlist(tapply(temp_MA$Variance,temp_MA$pop,function(x){
	sort(x)[c(0.05,0.95)*length(x)]})),nrow=2))


	
#	arrows(k, vect90[,1],k, vect90[,2],
#	code=3,length=.1,angle=90,col= "forestgreen")

	arrows(k+c(-.1,.1), vect90_MA[,1],k+c(-.1,.1), vect90_MA[,2],
	code=3,length=.1,angle=90,col= c("magenta","yellow"))


#points(post_modes_temp[8]~k,pch=16,col= "forestgreen")
points(post_modes_temp_MA~c(k+c(-.1,.1)),pch=16,col= c("magenta","yellow"))

}


legend(0,0.15,c("N2","PB306"),lwd=1.5,col=c("magenta","yellow"),bty="n",pch=16)


