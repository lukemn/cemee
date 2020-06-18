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

vect_col =c("black","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick","forestgreen")
vect_y_pops=rep(c(0,50,100,-33),c(1,3,3,1))

vect_margin=rbind(c(1.5,3,1.5,0),c(1.5,1.5,1.5,1.5),c(1.5,0,1.5,3),
c(2,3,1,0),c(2,1.5,1,1.5),c(2,0,1,3))
adjust_x=c(0,rep(c(-5,0,5),2),0)
#vect_margin=vect_margin/1.5

vect_main=c(expression(y[1]),expression(y[2]),expression(y[3]),
expression(y[4]),expression(y[5]),expression(y[6]))
par(mfrow=c(2:3))
for(k in 1:6){
par(mar= vect_margin[k,])
temp_pop=subset(list_pts, Rotated_traits==k & !is.na(Variance))
post_modes_temp=tapply(temp_pop$Variance,temp_pop$pop,posterior.mode)

plot(post_modes_temp[c(1:7,10)]~vect_y_pops,ylim=c(0.0008,0.7),type="n",log="y",xlim=c(-40,120),xlab="Generation",ylab="",bty="n",xaxt="n",yaxt="n",main= vect_main[k])

if(k>3) axis(side=1,at=c(0,50,100))
if(k %in% c(1,4)){
 axis(side=2,at=c(0.01,0.05,0.15,0.5))
}else{
 axis(side=2,at=c(0.01,0.05,0.15,0.5),labels=rep("",4))
}

vect90=t(matrix(unlist(tapply(temp_pop$Variance,temp_pop$pop,function(x){
	sort(x)[c(0.05,0.95)*length(x)]})),nrow=2))
	
	arrows(vect_y_pops+adjust_x, vect90[,1], vect_y_pops+adjust_x, vect90[,2],
	code=3,length=.1,angle=90,col= vect_col)

Ne_est=1000
lines(0:100,(1-2/(Ne_est))^(0:100)*median(post_modes_temp[1]),lwd=2,lty=2)
points(post_modes_temp[c(1:8)]~I(vect_y_pops+ adjust_x),pch=16,col= vect_col)



}





