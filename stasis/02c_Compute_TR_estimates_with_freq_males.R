rm(list=ls())
library(grImport)
library(boot)
library(data.table)
library(lme4)
library(plyr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(matrixStats)


## File that pre-process the population data (Males)

load("~/PATH/TO/DIR/Cemee_Pop_WI/State_frequencies.RData")
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")

#add the state frequencies to the final export DF
final_export=merge(final_export,state_final[,c("exper_data_id","data_id","data_group_name","mean_still","mean_forward","mean_backward")])

#### MD and WI ?
data_MD <- subset(final_export,data_group_name%in%"B310" & population2=="pops" &
                    tstrsplit(pop_label,"_")[[2]]=="male")
data_WI <- subset(final_export,population2=="WI")

##
data_populations0 <- subset(final_export,data_group_name%in%paste0("B30",1:9) & population2=="pops" & tstrsplit(pop_label,"_")[[2]]=="male")
data_populations=rbind(data_populations0, data_MD, data_WI)
dim(data_populations) #168 measurements

data_populations$population=tstrsplit(data_populations$population,"noM")[[1]]
data_populations$temperature = (data_populations$temperature-mean(data_populations$temperature))/sd(data_populations$temperature)
data_populations$rel_humidity = (data_populations$rel_humidity-mean(data_populations$rel_humidity))/sd(data_populations$rel_humidity)
data_populations$logD = (data_populations$logD-mean(data_populations$logD))/sd(data_populations$logD)
vect_P_traits2=paste0(vect_P_traits,"_pred_02")

## A linear model to estimate the mean value for each population
vect_traits_labels <- paste0(rep(c("Still","Forward","Backward"),each=2),'-',c("Forward","Backward","Still","Backward","Still","Forward"))
vect_traits_labels=c(vect_traits_labels,rep("Mean trait frequency",3))
vect_traits_labels_short=c(paste0(rep(c("S","F","B"),each=2),c("F","B","S","B","S","F")),c("Still","Forward","Backward"))

data_populations$population <- as.factor(as.character(data_populations$population))
beta=array(,c(9,49))
diff_G=array(,c(9,49))

CI_list=list()
vect_P_traits2=c(vect_P_traits2,"mean_still","mean_forward","mean_backward")

ylim_vect <- cbind(c(-2.2,-2.5,-2.4,-5.5,-1.3,-3.4,-2.6,-1,-4),
                   c(0,-.3,0,-2,.8,-.7,.5,2.5,-1))

k=0
for(i in 1:length(vect_P_traits2)){
  k=k+1
  k2=0
  
  temp_mod1 <- lmer(data_populations[, vect_P_traits2[i]]~population-1+(1| data_group_name),data= data_populations)
  
  temp_CI <- confint(temp_mod1)
  temp_CI <- temp_CI[3:nrow(temp_CI),]
  CI_list[[k]]=temp_CI
  mean_P <- as.numeric(summary(temp_mod1)$coef[,1])#[1:50]
  pop_lev= levels(data_populations$population)
  vA1 <- c("OF5","A0","A110","A130","A160","A1100");va1_gen=c(-5,0,10,30,60,100)
  
  temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
  temp_P=temp_P[order(temp_P[,2]),]
  
  for(k3 in 2:length(vA1)){
    k2=k2+1
    beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
    diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
  }
  
  for(temp_i in 2:3){
    va1_gen=c(0,10,30,60,100);vA1 = paste0("A", temp_i,va1_gen); vA1[1]="A0"
    temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
    temp_P=temp_P[order(temp_P[,2]),]
    
    for(k3 in 2:length(vA1)){
      k2=k2+1
      beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
      diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
    }
  }
  
  temp_i=4
  va1_gen=c(0,10,30,60,100,140);vA1 = paste0("A", temp_i,va1_gen); vA1[1]="A0"
  temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
  temp_P=temp_P[order(temp_P[,2]),]
  
  for(k3 in 2:length(vA1)){
    k2=k2+1
    beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
    diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
  }
  
  for(temp_i in 5:6){
    va1_gen=c(0,10,60,100,140);vA1 = paste0("A", temp_i,va1_gen); vA1[1]="A0"
    temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
    temp_P=temp_P[order(temp_P[,2]),]
    
    for(k3 in 2:length(vA1)){
      k2=k2+1
      beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
      diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
    }
  }
  
  ### CA populations
  for(temp_i in 1:3){
    va1_gen=c(0,5,10,36,50,68,100);vA1 = paste0("CA", temp_i,va1_gen)
    vA1[1]="A6140"; va1_gen= va1_gen+140
    
    temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
    temp_P=temp_P[order(temp_P[,2]),]
    
    for(k3 in 2:length(vA1)){
      k2=k2+1
      beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
      diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
    }
  }
  for(temp_i in 4:5){
    va1_gen=c(0,32,66);vA1 = paste0("CA", temp_i,va1_gen)
    vA1[1]="A6140"; va1_gen= va1_gen+140
    
    temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
    temp_P=temp_P[order(temp_P[,2]),]
    
    for(k3 in 2:length(vA1)){
      k2=k2+1
      beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
      diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
    }
  }
  va1_gen=c(140,172);vA1 = c("A6140","CA632")
  temp_P= cbind(mean_P[pop_lev%in%vA1], va1_gen[order(vA1)])
  temp_P=temp_P[order(temp_P[,2]),]
  
  for(k3 in 2:length(vA1)){
    k2=k2+1
    beta[k,k2]= temp_P[c(k3),1]-temp_P[c(k3-1),1]
    diff_G[k,k2]= temp_P[c(k3),2]-temp_P[c(k3-1),2]
    
  }
}

#######
for(j in 1:ncol(beta)){
  for(i in 1:nrow(beta)){
    beta[i,j]=beta[i,j]/diff_G[i,j]
  }}

beta_coordinates=rbind(
  c(1:4,1,6:8,1,10:12,1,14:17,1,19:21,1,23:25,26,27:31,26,33:37,26,39:43,26,45,26,47,26),
  c(2:49))

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq_Males.RData")
