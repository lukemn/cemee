rm(list = ls())
gc()
library(MCMCglmm)
library(psych)
library(ggplot2)
library(dplyr)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(Rmisc)
library(dae)
library(nlme)
library(parallel)
library(RColorBrewer)

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}

load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData')

vect_Pops=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")

####
MCMCtot <- nrow(VCV_mat[[1]]$VCV_Mat)
MCMCsamp <- 1000 
n <- 6 #number of traits
m <- 7 #number of matrices to compare
r <- 3 #number of random effects specified in the model.
traitnames <- vect_P_traits #trait names
Gnames <- vect_Pops

MCMCarray <- array(, c(MCMCsamp, (n^2) * r, m)) #empty array
MCMCarray[, , 1] <- as.matrix(VCV_mat[[1]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 2] <- as.matrix(VCV_mat[[2]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 3] <- as.matrix(VCV_mat[[3]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 4] <- as.matrix(VCV_mat[[4]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 5] <- as.matrix(VCV_mat[[5]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 6] <- as.matrix(VCV_mat[[6]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 7] <- as.matrix(VCV_mat[[7]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])

Garray <- array(, c(n, n, m, MCMCsamp))
dimnames(Garray) <- list(traitnames, traitnames, Gnames)
Parray <- array(, c(n, n, m, MCMCsamp))
dimnames(Parray) <- list(traitnames, traitnames, Gnames)

Earray1 <- array(, c(n, n, m, MCMCsamp))
dimnames(Earray1) <- list(traitnames, traitnames, Gnames)
Earray2 <- array(, c(n, n, m, MCMCsamp))
dimnames(Earray2) <- list(traitnames, traitnames, Gnames)

for (i in 1:m) {
	for (j in 1:MCMCsamp) {
		G <- matrix(MCMCarray[j, 1:(n^2), i], ncol = n)
		R1 <- matrix(MCMCarray[j, ((n^2) + 1):((n^2) * 2), i], ncol = n)
		R2 <- matrix(MCMCarray[j, (((n^2) * 2) + 1):((n^2) * 3), i], ncol = n)
		Garray[, , i, j] <- G
		Earray1[, , i, j] <- R1
		Earray2[, , i, j] <- R2	
		Parray[, , i, j] <- G + R1 + R2
	}
}

source('~/PATH/TO/DIR/Tensor_Analysis/functions_tensor.R', chdir = TRUE)

HHGarray <- array(, c(n, n, m, MCMCsamp))
for (k in 1:MCMCsamp) {
	for (j in 1:m) {
		P <- inv.rootP(Parray[, , j, k])
		HHGarray[, , j, k] <- P %*% Garray[, , j, k] %*% P
	}
}

## Create a pedigree

df_for_tensor= data_for_G_estimation

ped_all = rbind(
#first all the lines
data.frame(id = as.character(unique(df_for_tensor$pop_label)), dam = NA, sire = NA,stringsAsFactors=FALSE),
#then all the phenotyped lines
data.frame(id = 1:nrow(df_for_tensor), dam = as.character(df_for_tensor$pop_label), sire = as.character(df_for_tensor$pop_label),stringsAsFactors=FALSE)
)
for(i in 1:3) ped_all[,i]=as.factor(ped_all[,i])

population_for_ped <- data.frame(pop_label= c(as.character(unique(df_for_tensor$pop_label)),as.character(df_for_tensor$pop_label)),stringsAsFactors=FALSE)
population_for_ped$population=NA
for(i in 1:nrow(population_for_ped)) population_for_ped$population[i] = as.character(subset(unique(df_for_tensor[,c("population","pop_label")]),pop_label==population_for_ped$pop_label[i])$population)

rand.Garray <- array(, c(n, n, m, MCMCsamp))
rand.Garray_corrected <- array(, c(n, n, m, MCMCsamp))

dimnames(rand.Garray) <- list(traitnames, traitnames, Gnames)

df_for_tensor$population=as.factor(as.character(df_for_tensor$population))
rm(i)

# Here we save a file that could be used on a server to compute the randomized
# eigentensors that are computationaly demanding.
save(list=ls(),file='~/PATH/TO/DIR/Tensor_Analysis/File_for_parallel_processing.RData')

run_parallel_MCMC <- function(i){

	library(MCMCglmm)
	library(dae)
	load('~/PATH/TO/DIR/Tensor_Analysis/File_for_parallel_processing.RData')
	rand.Garray_corrected_parallel <- array(, c(n, n, m))
	
	A6140.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[1]), Garray[, , 1, i]/2)

	CA150.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[2]), Garray[, , 2, i]/2)
	CA250.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[3]), Garray[, , 3, i]/2)
	CA350.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[4]), Garray[, , 4, i]/2)

	
	CA1100.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[5]), Garray[, , 5, i]/2)
	CA2100.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[6]), Garray[, , 6, i]/2)
	CA3100.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[7]), Garray[, , 7, i]/2)


	a.pop <- cumsum(as.numeric(tapply(ped_all$id, population_for_ped$population, length))[c(1,3,5,7,2,4,6)])
	pop.bv <- rbind(A6140.bv, CA150.bv, CA250.bv, CA350.bv, CA1100.bv, CA2100.bv, CA3100.bv)

	rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1], replace = F), ]
	
	## Here we have to compute the random Garray using the morissey technique and save them in a list
	
	ped_lines_all <- subset(ped_all,!is.na(dam))
	rand.model_MCMC=list()
	for(k in 1:7){

	k_pop=Gnames[k]
	ped_lines_current <- subset(ped_lines_all, df_for_tensor$population==k_pop)
	if(k %in% 1:4) sire.bvs <- rand.pop.bv[substring(ped_all$id,1,5)==k_pop & is.na(ped_all$dam),]
	if(k %in% 5:7) sire.bvs <- rand.pop.bv[substring(ped_all$id,1,6)==k_pop & is.na(ped_all$dam),]
	# Add the proper sire label, that will differ from the row name (because of the shuffling)
	if(k %in% 1:4) sire.bvs <- cbind(data.frame(sire=ped_all$id[substring(ped_all$id,1,5)==k_pop & is.na(ped_all$dam)]),sire.bvs)
	if(k %in% 5:7) sire.bvs <- cbind(data.frame(sire=ped_all$id[substring(ped_all$id,1,6)==k_pop & is.na(ped_all$dam)]),sire.bvs)	
	
	ped_lines_current <- merge(ped_lines_current,sire.bvs)[,c(1,4:9)]
	#Vectors of E variance	
	z <- t(apply(ped_lines_current[,2:7],1,function(x){rmvnorm(x+rep(0,6),Earray2[,,k,i])}))
	ped_lines_current[,2:7] <- z
	names(ped_lines_current) <- c("pop_label",traitnames)
	phen.var = diag(nb_trait) * diag(var(z))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait)), R = list(V = phen.var/3, n = nb_trait))


	rand.model_MCMC.temp <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  trait - 1, random = ~us(trait):pop_label ,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = ped_lines_current, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000,thin=100)
	
		rand.Garray_corrected_parallel[,,k] <- matrix(posterior.mode(rand.model_MCMC.temp$VCV)[1:36], ncol = n)
	}
	
return(rand.Garray_corrected_parallel)	
}


clust <- makeCluster(25)
param_list=list()
for(i in 1:1000) param_list[[i]] <- i

List_output <-parLapply(clust, param_list , run_parallel_MCMC)

save(list=ls(),file="~/PATH/TO/DIR/Tensor_Analysis/Tensor_processed.Rdata")

### End of the parallel thread

#### Back on local computer
load('~/PATH/TO/DIR/Tensor_analysis/Tensor_processed.Rdata')
source('~/PATH/TO/DIR/Tensor_analysis/functions_tensor.R', chdir = TRUE)


for(i in 1:MCMCsamp){
	for(k in 1:m){
		 rand.Garray_corrected[,,k,i] <- matrix(List_output[[i]][,,k], ncol = n)
	}
}
dimnames(rand.Garray_corrected) <- list(traitnames, traitnames, Gnames)
MCMC.covtensor <- covtensor(Garray)
#MCMC.covtensor$tensor.summary
nnonzero <- min(n * (n + 1)/2, m - 1)
MCMC.covtensor.rand <- covtensor(rand.Garray_corrected)


HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.95), 
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.95))
round(HPD.eT.val, 3)


HPD.eT.val_80 <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.80), 
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.80))
round(HPD.eT.val_80, 3)

#
# Figure A1
quartz()
par(mfrow=c(1,1))

plot((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),xlab="",ylab=expression(alpha),pch=16,cex=1,xaxt="n",frame.plot=F,xlim=c(0.5,6.5),ylim=c(0,max(HPD.eT.val)),type="n")
axis(1,at=1:nnonzero,labels=c(paste("E",rep(1:nnonzero),sep="")))


arrows((1:nnonzero)-0.2, HPD.eT.val[,2],(1:nnonzero)-0.2,HPD.eT.val[,1],length=0.1,angle=90,code=3)
arrows((1:nnonzero)-0.2, HPD.eT.val_80[,2],(1:nnonzero)-0.2, HPD.eT.val_80[,1],length=0.1,angle=90,code=0,lwd=4,col="red")

arrows((1:nnonzero)+0.2, HPD.eT.val[,4],(1:nnonzero)+0.2,HPD.eT.val[,3],length=0.1,angle=90,lty=5,code=3)

arrows((1:nnonzero)+0.2, HPD.eT.val_80[,4],(1:nnonzero)+0.2, HPD.eT.val_80[,3],length=0.1,angle=90,code=0,lwd=4,col="orange")

points((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),pch=21,cex=1,bg="black")

points((1:nnonzero)+0.2, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),pch=8,cex=1,col="grey")

legend(3.5,.013,legend=c("observed","randomised"),lty=c(1,5),pch=c(16,8),cex=1,bty="n")


HPD.tensor.coord <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord) <- list(Gnames,c("lower","upper"), paste("E",1:6,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.95)[1:2]
  }
}


HPD.tensor.coord_80 <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord_80) <- list(Gnames,c("lower","upper"), paste("E",1:6,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord_80[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.80)[1:2]
  }
}





#Figure A2
par(mfrow=c(1,2))

for (k in 1:2){  
plot(1:m,MCMC.covtensor$av.G.coord[,k,],ylab="",xlab="",pch=16,xaxt="n",frame.plot=F,xlim=c(0.5,m+.5),ylim=c(-.3,.5),main = "",type="n")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.tensor.coord[,2,k],1:m,HPD.tensor.coord[,1,k],length=0.05,angle=90,code=3)

arrows(1:m,HPD.tensor.coord_80[,2,k],1:m, HPD.tensor.coord_80[,1,k],length=0.05,angle=90,code=0,col="red",lwd=3)

points(1:m,MCMC.covtensor$av.G.coord[,k,],pch=21,bg="grey")
mtext(dimnames(MCMC.covtensor$av.G.coord)[[2]][k],side=3,at=0,font=2)
}

round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)

# How much variation is explained ?

abs(MCMC.covtensor$tensor.summary[1,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))
#e11: 66%
abs(MCMC.covtensor$tensor.summary[2,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))
#e12:18%



abs(MCMC.covtensor$tensor.summary[7,2])/sum(abs(MCMC.covtensor$tensor.summary[(n+1):(2*n),2]))
#e11: 60%
abs(MCMC.covtensor$tensor.summary[8,2])/sum(abs(MCMC.covtensor$tensor.summary[(n+1):(2*n),2]))
#e12:31%




e11 <- c(as.numeric(MCMC.covtensor$tensor.summary[1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e12 <- c(as.numeric(MCMC.covtensor$tensor.summary[2,3:dim(MCMC.covtensor$tensor.summary)[2]]))

e13 <- c(as.numeric(MCMC.covtensor$tensor.summary[3,3:dim(MCMC.covtensor$tensor.summary)[2]]))

e11.proj <- apply(Garray, 3:4, proj, b = e11)
e12.proj <- apply(Garray, 3:4, proj, b = e12)
e13.proj <- apply(Garray, 3:4, proj, b = e13)
HPD.e11 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95)
HPD.e12 <- HPDinterval(t(as.mcmc(e12.proj)),prob = 0.95)
HPD.e13 <- HPDinterval(t(as.mcmc(e13.proj)),prob = 0.95)

e21 <- c(as.numeric(MCMC.covtensor$tensor.summary[n+1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e22 <- c(as.numeric(MCMC.covtensor$tensor.summary[n+2,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e23 <- c(as.numeric(MCMC.covtensor$tensor.summary[n+3,3:dim(MCMC.covtensor$tensor.summary)[2]]))

e21.proj <- apply(Garray, 3:4, proj, b = e21)
e22.proj <- apply(Garray, 3:4, proj, b = e22)
e23.proj <- apply(Garray, 3:4, proj, b = e23)
HPD.e21 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.95)
HPD.e22 <- HPDinterval(t(as.mcmc(e22.proj)),prob = 0.95)
HPD.e23 <- HPDinterval(t(as.mcmc(e23.proj)),prob = 0.95)


HPD.e11_80 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.80)
HPD.e12_80 <- HPDinterval(t(as.mcmc(e12.proj)),prob = 0.80)

HPD.e21_80 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.80)
HPD.e22_80 <- HPDinterval(t(as.mcmc(e22.proj)),prob = 0.80)



par(mfrow=c(2,2))
par(mar=c(5,4,4,2))
plot(1:m,rowMeans(e11.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e11))),las=2,type="n")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.e11[,2],1:m,HPD.e11[,1],length=0.1,angle=90,code=3)
arrows(1:m,HPD.e11_80[,2],1:m,HPD.e11_80[,1],length=0.1,angle=90,code=0,lwd=3,col="red")
points(1:m,rowMeans(e11.proj),pch=21,bg="grey")
mtext("e11 (66%)",side=3,at=0,font=2)

plot(1:m,rowMeans(e12.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e12))),las=2,type="n")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.e12[,2],1:m,HPD.e12[,1],length=0.1,angle=90,code=3)
arrows(1:m,HPD.e12_80[,2],1:m,HPD.e12_80[,1],length=0.1,angle=90,code=0,lwd=3,col="red")
points(1:m,rowMeans(e12.proj),pch=21,bg="grey")
mtext("e12 (18%)",side=3,at=0,font=2)

plot(1:m,rowMeans(e21.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e21))),type="n")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.e21[,2],1:m,HPD.e21[,1],length=0.1,angle=90,code=3)
arrows(1:m,HPD.e21_80[,2],1:m,HPD.e21_80[,1],length=0.1,angle=90,code=0,lwd=3,col="red")
points(1:m,rowMeans(e21.proj),pch=21,bg="grey")
mtext("e21 (60%)",side=3,at=0,font=2)

plot(1:m,rowMeans(e22.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e22))),type="n")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.e22[,2],1:m,HPD.e22[,1],length=0.1,angle=90,code=3)
arrows(1:m,HPD.e22_80[,2],1:m,HPD.e22_80[,1],length=0.1,angle=90,code=0,lwd=3,col="red")
points(1:m,rowMeans(e22.proj),pch=21,bg="grey")
mtext("e22 (31%)",side=3,at=0,font=2)

# Now I need to compute confidence intervals for the correlations

# Sample in the posterior of the gamma est. MCMC model
vect_sampling <- 1:1000
cor_dist <- array(dim=c(4,length(vect_sampling),6))

for(i in 1:length(vect_sampling)){
# Pick a Gamma matrix from the posterior distribution
temp_vect =(model_MCMC$Sol[vect_sampling[i],2:22])[c(7:21,1:6)]

rdm.gamma <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

cor_dist[1,i,] <- cor(e11,eigen(rdm.gamma)$vectors)
cor_dist[2,i,] <- cor(e12,eigen(rdm.gamma)$vectors)
cor_dist[3,i,] <- cor(e21,eigen(rdm.gamma)$vectors)
cor_dist[4,i,] <- cor(e22,eigen(rdm.gamma)$vectors)

}
for(k in 1:4){
for(i in 1:ncol(cor_dist[k,,])){
	cor_dist[k,,i] <- sort(abs(cor_dist[k,,i]))
	}
}

vect_col <- brewer.pal(n = 6, name = 'Dark2')

par(mfrow=c(2,2))
for(k in 1:4){
plot(density(cor_dist[k,,6]),xlim=c(0,1),las=1,bty="n",col=vect_col[6],main="",xlab="Correlation")
for(i in 1:5) lines(density(cor_dist[k,,i]),col=vect_col[i])
legend(.2,c(15,6,15,6)[k],expression(y[1], y[2], y[3], y[4], y[5], y[6]),col=vect_col,lwd=1,bty='n',ncol=2,cex=.8)
}

save(list=ls(),file='~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_tensor.RData')

