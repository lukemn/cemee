rm(list=ls())
  

#### Fig. S4 & S5
for(male_herm in c(1:2)){
if(male_herm==1){
	load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")

	ylim_vect <- cbind(c(-2.2,-2.7,-2.4,-5.5,-1.3,-3.4,-2.5,-1.5,-4),
c(-0.6,-1.2,0,-2.9,.8,-1,.5,2.5,-2))
}else{
	load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq_Males.RData")

ylim_vect <- cbind(c(-1.5,-2,-2.4,-3.5,-1.3,-2,NA,NA,NA),
c(0.5,0,-0.2,-1,.8,-.1,NA,NA,NA))
}


uniq_TP=NULL
quartz(height =5, width =10)
par(mfrow=c(2,3))

k=0
for(i in 1:6){
k=k+1
k2=0

temp_mod1 <- lmer(data_populations[, vect_P_traits2[i]]~population-1+(1| data_group_name),data= data_populations)

mean_P <- as.numeric(summary(temp_mod1)$coef[,1])
pop_lev= levels(data_populations$population)

vA1 <- c("OF5","A0","A110","A130","A160","A1100");va1_gen=c(-5,0,10,30,60,100)

temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
temp_P=temp_P[order(temp_P[,2]),]

vect_ylab=""
if(i %in%c(1,4)) vect_ylab="log transition rate"
plot(temp_P[,1]~temp_P[,2],type="b",xlim=c(-40,240),ylim=c(ylim_vect[i,]),main= paste0(vect_traits_labels[i]," (",vect_traits_labels_short[i],")"),xaxt="n",bty="n",xlab="Generations",ylab= vect_ylab,las=1,cex.lab=1.2)

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


arrows(0,CI_list[[k]][1,1],0,CI_list[[k]][1,2],code=3,angle=90,length=.05,lwd=2)
arrows(140,CI_list[[k]][25,1],140,CI_list[[k]][25,2],code=3,angle=90,length=.05,lwd=2)
points(temp_P[,1]~temp_P[,2],type="p",cex=1.3,pch=21,bg="white",lwd=2)

axis(side=1,at=c(-33,-5,0,10,30,60,100,140,145,150,172,176,190,206,208,240),labels=c("-33","","0","","","60","","140 [0]",rep("",4),"190 [50]",rep("",2),"240 [100]"))

}
}