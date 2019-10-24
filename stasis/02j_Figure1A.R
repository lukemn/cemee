rm(list=ls())
library(data.table)
library(lme4)
library(plyr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(matrixStats)
library(boot) 

ylim_vect_2=rbind(c(0,.65),c(.25,1),c(0,.1))


for(male_herm in c(1:2)){

if(male_herm==1) load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")

if(male_herm==2) load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq_Males.RData")

i=7;k=i;k2=0
vect_ylab="Still frequency"
if( male_herm==1) par(mar=c(5,6,2,3))


if(male_herm==1) plot(1,1,type="n",xlim=c(-60,240),ylim= ylim_vect_2[1,],main= "",xaxt="n",bty="n",xlab="",ylab= vect_ylab,las=1,yaxt="s",cex.lab=1.3)

if(male_herm==2) points(1,1,type="n",xlim=c(-60,240),ylim= ylim_vect_2[1,],main= "",xaxt="n",bty="n",xlab="",ylab= "",las=1,yaxt="n",cex.lab=1.3)
if(male_herm==1){
ylimit_polygon= rep(ylim_vect_2[1,],each=2)
polygon(x=c(-63,-3,-3,-63),y= ylimit_polygon,col=adjustcolor("darkgreen", alpha.f = 0.3),border=NA)
polygon(x=c(143,245,245,143),y= ylimit_polygon,col=adjustcolor("dodgerblue2", alpha.f = 0.3),border=NA)
}
	vect_phenotypes= data_populations[, "mean_still"]
	temp_mod1 <- lmer(vect_phenotypes ~population-1+(1| data_group_name),data= data_populations)
	mean_P <- as.numeric(summary(temp_mod1)$coef[,1])
	mean_P=inv.logit(mean_P)


pop_lev= levels(data_populations$population)


vA1 <- c("OF5","A0","A110","A130","A160","A1100");va1_gen=c(-5,0,10,30,60,100)

temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]

points(temp_P[,1]~temp_P[,2],type="b")

if(male_herm==1){
v_WI=pop_lev%in%vect_WI
for(i_WI in c(1:length(v_WI))[v_WI]){
points(c(mean_P[i_WI],mean_P[pop_lev =="OF5"])~c(-33,-5),pch=21,bg=c("pink","white"),type="b")
}
}
for(temp_i in 2:3){
va1_gen=c(0,10,30,60,100);vA1 = paste0("A", temp_i,va1_gen); vA1[1]="A0"
temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="b")
}

temp_i=4
va1_gen=c(0,10,30,60,100,140);vA1 = paste0("A", temp_i,va1_gen); vA1[1]="A0"
temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="b")

for(temp_i in 5:6){
va1_gen=c(0,10,60,100,140);vA1 = paste0("A", temp_i,va1_gen); vA1[1]="A0"
temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="b")
}


### CA populations
for(temp_i in 1:3){
va1_gen=c(0,5,10,36,50,68,100);vA1 = paste0("CA", temp_i,va1_gen)
vA1[1]="A6140"; va1_gen= va1_gen+140

temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="b")
}

for(temp_i in 4:5){
va1_gen=c(0,32,66);vA1 = paste0("CA", temp_i,va1_gen)
vA1[1]="A6140"; va1_gen= va1_gen+140

temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="b")
}
va1_gen=c(140,172);vA1 = c("A6140","CA632")
temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="b")


## bigger points for A0 and A6140
va1_gen=c(0,140);vA1 = c("A0","A6140")
temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]


arrows(0,inv.logit(CI_list[[k]])[1,1],0,inv.logit(CI_list[[k]])[1,2],code=3,angle=90,length=.05,lwd=2)
arrows(140,inv.logit(CI_list[[k]])[25,1],140,inv.logit(CI_list[[k]])[25,2],code=3,angle=90,length=.05,lwd=2)
points(temp_P[,1]~temp_P[,2],type="p",cex=1.3,pch=21,bg=c("white",grey(.5)),lwd=2)

## CA[1-3]50 and CA[1-3]100
va1_gen=rep(c(190,240),each=3);vA1 = c("CA150","CA250","CA350","CA1100","CA2100","CA3100")
temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]
points(temp_P[,1]~temp_P[,2],type="p",pch=21,bg=c("cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick","forestgreen"))

}
mtext(side=1,"Generations",padj=3.5,cex=1.2)	



axis(side=1,at=c(-33,0,60,140,190,240),labels=c("-33","0","60","140 [0]","190 [50]","240 [100]"),cex.axis=.9,line=NA,tick=FALSE)
axis(side=1,at=c(-33,-5,0,10,30,60,100,140,145,150,172,176,190,206,208,240),labels=NA)


points(90,.05,pch=-0x2642L,cex=2.5)

points(90,.55,pch=-0x2642L,cex=2.5)
points(88.53,.5195,pch=-0x2640L,cex=2.5,col="black")





