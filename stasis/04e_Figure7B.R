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

# Rotate all the population and lines to prepare for further analysis
names(mean_phen_values)[1]<-'pop_label'
mean_phen_values=merge(mean_phen_values,unique(final_export[,c("pop_label","population",'population2')]))


A= VCV_mat[[1]]$G1_mat
proj_M <- eigen(gamma)$vectors
## Rotate the points

rotated_phenotypes <- NULL
for(i in 1:nrow(final_fertility)) rotated_phenotypes <- rbind(rotated_phenotypes,t(t(proj_M)%*%t(final_fertility[i,9:14])))

rotated_phenotypes_all_A6 <- subset(mean_phen_values,substring(pop_label,1,2)=="A6" & !population2%in%c("pops","WI"))[,2:7]

for(i in 1:6) rotated_phenotypes_all_A6[,i]= rotated_phenotypes_all_A6[,i] - sav_Means_P_values[i]
for(i in 1:nrow(rotated_phenotypes_all_A6)) rotated_phenotypes_all_A6[i,] <- t(t(proj_M)%*%t(rotated_phenotypes_all_A6[i,]))


rotated_phenotypes_all_CA <- subset(mean_phen_values,substring(pop_label,1,2)=="CA" & !population2%in%c("pops","WI"))[,2:7]


for(i in 1:6) rotated_phenotypes_all_CA[,i]= rotated_phenotypes_all_CA[,i] - sav_Means_P_values[i]
for(i in 1:nrow(rotated_phenotypes_all_CA)) rotated_phenotypes_all_CA[i,] <- t(t(proj_M)%*%t(rotated_phenotypes_all_CA[i,]))


pops_herm=subset(mean_phen_values,population2=="pops" & tstrsplit(pop_label,"_")[[2]]=="herm" &
is.na(tstrsplit(pop_label,"noM")[[2]]) & !substring(pop_label,1,2)%in%c("CD","CM"))
sav_Means_P_values=as.numeric(sav_Means_P_values)
## Retain A0, A6140 and all CA
population_means_centered_for_plot <- t(pops_herm[,2:7])
pops_retained=pops_herm$population
for(i in 1:ncol(population_means_centered_for_plot)){
	population_means_centered_for_plot[,i] <- population_means_centered_for_plot[,i] -sav_Means_P_values
} 
colnames(population_means_centered_for_plot) <- pops_retained

### Wild isolates
WI_centered_for_plot <- t(subset(mean_phen_values,population2=="WI")[,2:7])
WI_retained=tstrsplit(subset(mean_phen_values,population2=="WI")$pop_label,"L0")[[1]]
for(i in 1:ncol(WI_centered_for_plot)){
	WI_centered_for_plot[,i] <- WI_centered_for_plot[,i] -sav_Means_P_values
} 
colnames(WI_centered_for_plot) <- WI_retained


## Wild Isolates

rotated_WI_means <- NULL

for(i in 1:ncol(WI_centered_for_plot)){
rotated_WI_means <- cbind(rotated_WI_means,as.numeric(t(t(proj_M)%*%(WI_centered_for_plot[,i]))))
}
colnames(rotated_WI_means)<-colnames(WI_centered_for_plot)



#### And now the figure 7B

focal_traits <- c(1,6)

scale_fert <- log(final_fertility$fertility)
scale_fert <- (scale_fert - mean(scale_fert))
scale_fert <- scale_fert+abs(min(scale_fert))+1
scale_fert <- 256-(round(256*(scale_fert-min(scale_fert))/max(scale_fert)))

library(gplots)
pts_phenotyped=cbind(rotated_phenotypes[, focal_traits[1]],rotated_phenotypes[,focal_traits[2]])
pts_mM<-rbind(colMins(pts_phenotyped),colMaxs(pts_phenotyped))
mult_fact=1.1

pts_mM=pts_mM*mult_fact
pts_mM[2,2]=.35
pts_mM[1,1]=-.42
fit2=NULL;ndim_z=50
Z1=seq(pts_mM[1,1], pts_mM[2,1],length.out= ndim_z)
Z2=seq(pts_mM[1,2], pts_mM[2,2],length.out= ndim_z)


compute_fit2 <- function(x){
	x=as.matrix(x)
	1/2*t(x)%*%diag(eigen(gamma)$values)%*%x
	}

for(i in Z1){
ZF=cbind(rep(i, ndim_z),rep(0,ndim_z),rep(0,ndim_z),rep(0,ndim_z),rep(0,ndim_z),Z2)
fit2 = rbind(fit2,apply(ZF, 1, compute_fit2 ))
}

vect_scale_color=max(abs(fit2))




plot_all_data=function(){

axis(side=1,at=c(-.5,0,.5));axis(side=2,at=c(-.5,0,.5))

abline(h=0,lty=2);abline(v=0,lty=2)

points(rotated_phenotypes[, focal_traits[1]],rotated_phenotypes[,focal_traits[2]],pch=21,cex=.8,bg=gray(seq(.4,.9,.4/256))[scale_fert],asp=1,col=grey(.3))


for(i in 1:ncol(population_means_centered_for_plot)){
rotated_pop_herm <- as.numeric(t(t(proj_M)%*%(population_means_centered_for_plot[,i])))
points(rotated_pop_herm[focal_traits[2]]~ rotated_pop_herm[focal_traits[1]],bg="yellow",pch=22)
}
}

vect_levels=seq(min(fit2),max(fit2),length.out=20)
vect_levels
vect_cols=bluered(32)[c(1:16,19,21,22)]

filled.contour(Z1,Z2,fit2,col= vect_cols,levels=vect_levels,xlim= pts_mM[,1],ylim= pts_mM[,2],plot.axes={plot_all_data()},xlab=expression(y[1]),ylab=expression(y[6]))






