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

#Produce Fig S12

final_merged =read.table("~/PATH/TO/DIR/Cemee_Pop_WI/Final_merged_data_NGM.txt",h=TRUE,sep="\t")
state_freq=read.table("~/PATH/TO/DIR/Cemee_Pop_WI/State_freq_NGM.txt",sep='\t',h=TRUE)

A6140_all=subset(final_merged,pop_label=="A6140L0_herm")
A6140_all_state=subset(state_freq,pop_label=="A6140L0_herm")
final_merged=subset(final_merged,data_group_name%in%c(paste0("B30",1:9),"B400"))

#### 

final_merged$WI_ID=tstrsplit(final_merged$data_id,"L")[[1]]
final_merged$WI_ID[final_merged$WI_ID==""]="LSJ2"
vect_WI_withB400 <- c("AB1","CB4507", "CB4852","CB4855" , "CB4856" ,"CB4858", "JU319","JU345" , "JU400" , "MY1","MY16","N2anc","PB306" ,"PX174","PX179" , "RC301","LSJ2","CX12311")
final_merged=subset(final_merged,WI_ID%in% vect_WI_withB400)

final_all=merge(final_merged,state_freq)
head(final_all[,c(16:21,23:25)])
vect_y_lim=c(min(final_all[,c(16:21,23:25)]),max(final_all[,c(16:21,23:25)]))


N2=subset(final_all, WI_ID=="N2anc")
CX12311=subset(final_all, WI_ID=="CX12311")
CB4507 =subset(final_all, WI_ID=="CB4507")

## WI - non-N2=
WInonN2=subset(final_all, data_group_name!="B400" & !WI_ID%in%c("N2anc","CB4507"))

traitM_nonN2=as.numeric(colMeans(WInonN2[,23:25]))
traitSd_nonN2=sqrt(colVars(as.matrix(WInonN2[,23:25])))

### Logit transformed

A6140_all=merge(A6140_all, A6140_all_state)

par(mfrow=c(1,3))
for(i in 1:3){

plot((traitM_nonN2)[i]~c(0),pch=21,ylim=cbind(c(-3.5,2),c(-2,3),c(-4,-2))[,i],bty="n",ylab="",xaxt="n",xlab="",cex=2,cex.lab=1.4,xlim=c(-.1,.12),main=c("still","forward","backward")[i],lwd=2)

if(i==1){
	mtext(side=2,"Frequency (logit)",padj=-3)
	}

	axis(side=1,at=c(-.025,.1),labels=c('Wild \nisolates',"A6140"),las=2,cex.axis=1.3)

arrows(0,traitM_nonN2[i]-traitSd_nonN2[i],0,traitM_nonN2[i]+traitSd_nonN2[i],lwd=1.5,code=3,angle=90,length=.07)

points((N2[,i+22])~ rep(-.05,nrow(N2)),bg="green",pch=21,cex=2)
points((CB4507[,i+22])~rep(-.05,nrow(CB4507)),bg="yellow",pch=21,cex=2)
points((CX12311[,i+22])~rep(-.05,nrow(CX12311)),bg="red",pch=21,cex=2)

points(mean(A6140_all[,i+21])~c(.1),col="black",pch=8,cex=2,lwd=1.2)
points(mean(A6140_all[,i+21])~c(.1),col="black",pch=1,cex=2.5,lwd=1.2)

points(WInonN2[,i+22]~jitter(c(rep(0,nrow(WInonN2)))),cex=1,col="grey",pch=21)

if(i==3) legend(0,-2,c('N2','CX12311'),pt.bg=c("green","red"),pch=21,pt.cex=1.2)

}
