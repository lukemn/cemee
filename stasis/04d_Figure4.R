rm(list = ls())
library(data.table)
library(pracma)
library(matrixStats)
library(ggplot2)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")
load("~/PATH/TO/DIR/Simulations/Main_model_output_list_with_metadata.RData")

#Layout from old figure
lmat=rbind(
c(0,0,6,7,0,rep(18,4)),
c(0,1,8,13,0,rep(18,4)),
c(0,2,9,14,0,rep(18,4)),
c(0,3,10,15,0,rep(18,4)),
c(0,4,11,16,0,0,19:21),
c(0,5,12,17,0,0,19:21))
lmat2=lmat[,6:9]-17
lmat2[lmat2<0]=0
lmat2=cbind(rep(0,6), lmat2)

#### The barplot now
eigen_values <- data.frame(pop=rep(vect_populations,each=6),gmax=c(
eigen(VCV_mat[[1]]$G1_mat)$values,
eigen(VCV_mat[[2]]$G1_mat)$values,
eigen(VCV_mat[[3]]$G1_mat)$values,
eigen(VCV_mat[[4]]$G1_mat)$values,
eigen(VCV_mat[[5]]$G1_mat)$values,
eigen(VCV_mat[[6]]$G1_mat)$values,
eigen(VCV_mat[[7]]$G1_mat)$values),numeig = rep(1:6,7),eigenV_norm=NA,pop2=rep(c(0,50,100),c(6,18,18)))
for(i in 1:7) eigen_values[eigen_values$pop== vect_populations[i],]$eigenV_norm=eigen_values[eigen_values$pop== vect_populations[i],]$gmax/sum(eigen_values[eigen_values$pop== vect_populations[i],]$gmax)

########################################################################

rdm_spl=sample(1: 10000)
variance_mat <- array(0,c(6,100,1))

for(i in 1:100){
		j=1
variance_mat[,i,j]=eigen(matrix(VCV_mat[[j]]$VCV_Mat[rdm_spl[i],1:36],6,6))$values
}

eigen_Sds=data.frame(pop="A6140", numeig=1:6,Sds=sqrt(rowVars(variance_mat[,,1])))
eigen_values2=merge(subset(eigen_values,pop=="A6140"), eigen_Sds)


for(k in c(5,1)){
tmp_sub=subset(eigen_values,substring(pop,4,4)==k)
tmp_CA=tapply(tmp_sub$eigenV_norm,tmp_sub$numeig,mean)
tmp_CA_sd=tapply(tmp_sub$eigenV_norm,tmp_sub$numeig,sd)

if(k==5){ eigen_values2=rbind(eigen_values2,data.frame(pop="CA1[1-3]50",gmax=NA, eigenV_norm=as.numeric(tmp_CA),Sds=as.numeric(tmp_CA_sd),pop2=50, numeig=1:6))
	}else{
eigen_values2=rbind(eigen_values2,data.frame(pop="CA1[1-3]100",gmax=NA, eigenV_norm=as.numeric(tmp_CA),Sds=as.numeric(tmp_CA_sd),pop2=100,numeig=1:6))
}
}


## We need to add the results of the simulation
simul_eigen=NULL
for(i in 1:length(G_mat_final)) simul_eigen=rbind(simul_eigen,eigen(G_mat_final[[i]])$values)
gmax_simul= colMedians(simul_eigen)
for(i in 1:nrow(simul_eigen)) simul_eigen[i,]=simul_eigen[i,]/sum(simul_eigen[i,])

eigen_values2=rbind(eigen_values2,
data.frame(pop='Genetic drift',numeig=1:6,gmax= gmax_simul, eigenV_norm =colMedians(simul_eigen),pop2=200,Sds=colSds(simul_eigen)))

layout(lmat2,widths=c(1,1,rep(4,3)), heights=c(.8,rep(2,5)))
par(mar=c(3,4,0,0),xaxt="s",yaxt="s")
vect_spaces=c(0,rep(rep(c(0,1.5),c(3,1)),6))[1:24]
barplot(eigen_values2[order(eigen_values2$numeig,eigen_values2$pop2),]$eigenV_norm,space= vect_spaces,
col=c("grey","dodgerblue2","firebrick3","green"),ylim=c(0,.75),yaxt="n")

axis(side=2,at=seq(0,.6,length.out=7),las=1,cex.axis=1.5)
mtext(side=3,"A",adj=-.09,padj=1.7,cex=2.1)
mtext(side=1,"Principal axes of genetic variance",padj=2,cex=1)
axis(side=1,line=NA,tick=FALSE,at=c(2,7.5,13,18.5,24,29.5),labels=c(expression(g[1]),expression(g[2]),expression(g[3]),expression(g[4]),expression(g[5]),expression(g[6])),cex.axis=1.4,padj=-.7)

mtext(side=2,"Proportion of genetic variance",padj=-4,adj=.3)

# Arrows for A6140
vect_arrows= cbind(subset(eigen_values2,pop=="A6140")$eigenV_norm-subset(eigen_values2,pop=="A6140")$Sds,
subset(eigen_values2,pop=="A6140")$eigenV_norm+subset(eigen_values2,pop=="A6140")$Sds)
vect_arrows_x=c(0.5,6,11.5,17,22.5,28)
arrows(vect_arrows_x,vect_arrows[,1], vect_arrows_x,vect_arrows[,2],code=3,angle=90,length=.05,lwd=1,col=grey(.5))

# Arrows for Drift
vect_arrows= cbind(subset(eigen_values2,pop=="Genetic drift")$eigenV_norm-subset(eigen_values2,pop=="Genetic drift")$Sds,
subset(eigen_values2,pop=="Genetic drift")$eigenV_norm+subset(eigen_values2,pop=="Genetic drift")$Sds)
vect_arrows_x=c(0.5,6,11.5,17,22.5,28)+3
arrows(vect_arrows_x,vect_arrows[,1], vect_arrows_x,vect_arrows[,2],code=3,angle=90,length=.05,lwd=1,col=grey(.5))

# Arrows CA 50 and CA 100
k=0
for(CA_nb in c(50,100)){
	k=k+1
vect_arrows= cbind(subset(eigen_values2,pop2== CA_nb)$eigenV_norm-subset(eigen_values2,pop2==CA_nb)$Sds,
subset(eigen_values2,pop2==CA_nb)$eigenV_norm+subset(eigen_values2,pop2==CA_nb)$Sds)
vect_arrows_x=c(0.5,6,11.5,17,22.5,28)+k
arrows(vect_arrows_x,vect_arrows[,1], vect_arrows_x,vect_arrows[,2],code=3,angle=90,length=.05,lwd=1,col="black")
}
legend(18,.58,bty="n",c("Ancestor (A6140)","Generation 50 (CA[1-3]50)","Generation 100 (CA[1-3]100)","Genetic drift (simulations)"),col=c(grey(.5),"dodgerblue2","firebrick3","green"),lwd=10,cex=1.4)



#### Angle for the first 3 EV
i1=1
par(yaxt="s",mar=c(0.2,0.7,0,0.2),xaxt="n")
for(i1 in 1:3){
A = VCV_mat[[1]]$G1_mat
eigVal_A <- eigen(A)$values; eigVec_A <- eigen(A)$vectors
yvect_gmax=cbind()

plot(c(-1, 1), c(0, 0), type = "b", lwd = 2, ylim = c(-1,1), asp = 1, 
		ylab = "", yaxt = "n",bty="n",las=1,xlab="")
if(i1==1){
	mtext(side=3,"B",adj=-.37,cex=2,padj=.5)
	mtext(side=2,expression(atop(paste("Orientation (",theta,") and"))))
	mtext(side=2,expression(paste("\n magnitude of ",g[i])),padj=-.6)
} 

axis(side=2,at=c(-1,0,1),labels=c("","",""),pos=-1.1)
for(j in 2:7){
temp_G <- VCV_mat[[j]]$G1_mat
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])

temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))
		points(c(-temp_coordinates[1], temp_coordinates[1]), c(-temp_coordinates[2], temp_coordinates[2]), type = "l", col = rep(c("dodgerblue2","firebrick3"),c(4,3))[j], lwd = 2)

}

for(temp_c in 1:length(G_mat_final)){

temp_G <- G_mat_final[[temp_c]]
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])
temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))
		points(c(-temp_coordinates[1], temp_coordinates[1]), c(-temp_coordinates[2], temp_coordinates[2]), type = "l", col = "green", lwd = 2)
}

if(i1==1) mtext(side=1,expression(paste(g[1]," (62%) ")),padj=-2,cex.axis=1.8)
if(i1==2) mtext(side=1,expression(paste(g[2]," (23%) ")),padj=-2,cex.axis=1.8)
if(i1==3) mtext(side=1,expression(paste(g[3]," (8%) ")),padj=-2,cex.axis=1.8)
}



