rm(list=ls())
library(data.table)
library(lme4)
library(plyr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(matrixStats)
library(boot) 
library(MCMCglmm)


load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")


### A6140
load('~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_A6140_.RData')

df_G1 <- NULL
for(i in 1:length(List_output)){
  df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
}

rm(List_output);gc()
save(list=ls(),file='~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_A6140_LIGHT.RData')

load('~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_A6140_LIGHT.RData')

plot_G_coefs <- function(true_G, df_G1,pop_name){
  par(mar=c(5,7,4,2))
  vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
  vProb <- .95
  
  
  vectX <- c(true_G$G1_mat)[vect_Var]
  plot(vectX,c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.16,.31),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
  mtext(side=2,"Phenotypic traits    \n Diagonal                                  Off-diagonal            ",padj=-2,cex=1.2)
  lines(c(0,0),c(24.5,8.5))
  lines(c(0,0),c(.5,5.5),col="red")
  axis(side=1,pos=0)
  axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                       "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                       "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                       "SF","SB","FS","FB","BS","BF"),las=1)
  
  temp_95 <- HPDinterval(true_G$VCV_Mat[,1:36],prob=.95)
  arrows(temp_95[vect_Var,1],c(24:10,6:1),temp_95[vect_Var,2],c(24:10,6:1),code=3,length=.05,angle=90)
  
  temp_80 <- HPDinterval(true_G$VCV_Mat[,1:36],prob=.8)
  arrows(temp_80[vect_Var,1],c(24:10,6:1), temp_80[vect_Var,2],c(24:10,6:1),code=3,length=0,angle=90,lwd=2,col="red")
  
  points(vectX,c(24:10,6:1),pch=21,bg="black",cex=.8)
  
  temp_95_rand <- HPDinterval(as.mcmc(df_G1[,1:36]),prob=.95)
  arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)-.3, temp_95_rand[vect_Var,2],c(24:10,6:1)-.3,code=3,length=.05,angle=90,col='grey')
  
  temp_80_rand <- HPDinterval(as.mcmc(df_G1[,1:36]),prob=.8)
  arrows(temp_80_rand[vect_Var,1],c(24:10,6:1)-.3, temp_80_rand[vect_Var,2],c(24:10,6:1)-.3,code=3,length=0,angle=90,lwd=2,col="orange")
  
  vectX_rand <- posterior.mode(as.mcmc(df_G1[,1:36]))[vect_Var]
  points(vectX_rand,c(24:10,6:1)-.3,pch=8,col="grey",cex=.8)
  
  legend(.1,23,c(paste0(pop_name," estimates"),"Null distribution"),pch=c(16,8),lwd=1)
}

plot_G_coefs(VCV_mat[[1]], df_G1,"A6140")

### CA- pops
for(popname in vect_populations[2:7]){
  
  load(paste0("~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_",popname,"_.RData"))
  
  df_G1 <- NULL
  for(i in 1:length(List_output)){
    df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
  }
  
  rm(List_output);gc()
  save(list=ls(),file=paste0("~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_",popname,"_LIGHT.RData"))
  
  rm(list=ls());gc
}

popname="A6140"
load(paste0("~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_",popname,"_LIGHT.RData"))

for(kpopname in 2:7){
  popname <- vect_populations[kpopname]
  load(paste0("~/PATH/TO/DIR/Random_G_Analysis_Cemee_Pop_WI_",popname,"_LIGHT.RData"))
  
  pdf(paste0("~/PATH/TO/DIR/Figures/",popname,"_matrix.pdf"))
  plot_G_coefs(VCV_mat[[kpopname]], df_G1,popname)
  
  dev.off()
}

### A6140 down-sampled vs all 

### A6140
load("~/PATH/TO/DIR/Processed_A6140_subset.RData")


vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95
true_G <- VCV_mat[[1]]
vectX <- c(true_G$G1_mat)[vect_Var]
plot(vectX,c(24:10,6:1),yaxt="n",bty="n",ylab="",xlim=c(-.16,.31),xlab="",xaxt="n",type='n')
lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")
axis(side=1,pos=0)
axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                     "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                     "SF","SB","FS","FB","BS","BF"),las=1)

temp_95 <- HPDinterval(true_G$VCV_Mat[,1:36],prob=.95)
arrows(temp_95[vect_Var,1],c(24:10,6:1),temp_95[vect_Var,2],c(24:10,6:1),code=3,length=.05,angle=90)

temp_80 <- HPDinterval(true_G$VCV_Mat[,1:36],prob=.8)
arrows(temp_80[vect_Var,1],c(24:10,6:1), temp_80[vect_Var,2],c(24:10,6:1),code=3,length=0,angle=90,lwd=2,col="red")

points(vectX,c(24:10,6:1),pch=21,bg="black",cex=.8)

temp_95_rand <- HPDinterval(as.mcmc(df_G2[,1:36]),prob=.95)
arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)-.3, temp_95_rand[vect_Var,2],c(24:10,6:1)-.3,code=3,length=.05,angle=90,col='grey')

temp_80_rand <- HPDinterval(as.mcmc(df_G2[,1:36]),prob=.8)
arrows(temp_80_rand[vect_Var,1],c(24:10,6:1)-.3, temp_80_rand[vect_Var,2],c(24:10,6:1)-.3,code=3,length=0,angle=90,lwd=2,col="orange")

vectX_rand <- posterior.mode(as.mcmc(df_G2[,1:36]))[vect_Var]
points(vectX_rand,c(24:10,6:1)-.3,pch=8,col="grey",cex=.8)

legend(.1,23,c(paste0(pop_name," estimates"),"Null distribution"),pch=c(16,8),lwd=1)

legend(.1,23,c(paste0("A6140 (189 RILs)"),paste0("A6140 (50 RILs)")),pch=c(16,8),lwd=1)


