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
library(corrplot)
library(Rmisc)
library(nlme)

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}



load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

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
for (i in 1:m) {
	for (j in 1:MCMCsamp) {
		G <- matrix(MCMCarray[j, 1:(n^2), i], ncol = n)
		R1 <- matrix(MCMCarray[j, ((n^2) + 1):((n^2) * 2), i], ncol = n)
		R2 <- matrix(MCMCarray[j, (((n^2) * 2) + 1):((n^2) * 3), i], ncol = n)
		Garray[, , i, j] <- G
		Parray[, , i, j] <- G + R1 + R2
	}
}

source('~/PATH/TO/DIR/functions_tensor.R', chdir = TRUE)

HHGarray <- array(, c(n, n, m, MCMCsamp))
for (k in 1:MCMCsamp) {
	for (j in 1:m) {
		P <- inv.rootP(Parray[, , j, k])
		HHGarray[, , j, k] <- P %*% Garray[, , j, k] %*% P
	}
}

## Create a pedigree

df_for_tensor= data_for_G_estimation

ped = data.frame(animal = df_for_tensor$pop_label, dam = NA, sire = NA)
rand.Garray <- array(, c(n, n, m, MCMCsamp))
dimnames(rand.Garray) <- list(traitnames, traitnames, Gnames)
df_for_tensor$population=as.factor(as.character(df_for_tensor$population))
for (i in 1:MCMCsamp) {

	A6140.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[1]), Garray[, , 1, i])
	CA150.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[2]), Garray[, , 2, i])
	CA250.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[3]), Garray[, , 3, i])
	CA350.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[4]), Garray[, , 4, i])

	CA1100.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[5]), Garray[, , 5, i])
	CA2100.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[6]), Garray[, , 6, i])
	CA3100.bv <- rbv(subset(ped, df_for_tensor$population == Gnames[7]), Garray[, , 7, i])


	a.pop <- cumsum(as.numeric(tapply(ped$animal, df_for_tensor$population, length))[c(1,3,5,7,2,4,6)])
	pop.bv <- rbind(A6140.bv, CA150.bv, CA250.bv, CA350.bv, CA1100.bv, CA2100.bv, CA3100.bv)
	rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1], replace = F), ]
	rand.Garray[, , 1, i] <- cov(rand.pop.bv[1:a.pop[1], ])
	rand.Garray[, , 2, i] <- cov(rand.pop.bv[(a.pop[1] + 1):a.pop[2], ])
	rand.Garray[, , 3, i] <- cov(rand.pop.bv[(a.pop[2] + 1):a.pop[3], ])
	rand.Garray[, , 4, i] <- cov(rand.pop.bv[(a.pop[3] + 1):a.pop[4], ])
	rand.Garray[, , 5, i] <- cov(rand.pop.bv[(a.pop[4] + 1):a.pop[5], ])
	rand.Garray[, , 6, i] <- cov(rand.pop.bv[(a.pop[5] + 1):a.pop[6], ])
	rand.Garray[, , 7, i] <- cov(rand.pop.bv[(a.pop[6] + 1):a.pop[7], ])
}
MCMC.covtensor <- covtensor(Garray)
#MCMC.covtensor$tensor.summary

nnonzero <- min(n * (n + 1)/2, m - 1)
MCMC.covtensor.rand <- covtensor(rand.Garray)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.9), 
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.9))

round(HPD.eT.val, 3)
#
# Figure A1
quartz()
par(mfrow=c(1,1))

plot((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),xlab="",ylab=expression(alpha),pch=16,cex=1,xaxt="n",frame.plot=F,xlim=c(0.5,7.5),ylim=c(0,max(HPD.eT.val)))
axis(1,at=1:nnonzero,labels=c(paste("E",rep(1:nnonzero),sep="")))
points((1:nnonzero)+0.2, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),pch=1,cex=1)
arrows((1:nnonzero)-0.2, unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)-0.2,HPD.eT.val[,1],length=0.1,angle=90)
arrows((1:nnonzero)-0.2, unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)-0.2,HPD.eT.val[,2],length=0.1,angle=90)
arrows((1:nnonzero)+0.2, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)+0.2,HPD.eT.val[,3],length=0.1,angle=90,lty=5)
arrows((1:nnonzero)+0.2, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)+0.2,HPD.eT.val[,4],length=0.1,angle=90,lty=5)
legend(3.5,.013,legend=c("observed","randomised"),lty=c(1,5),pch=c(16,1),cex=1,bty="n")


HPD.tensor.coord <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord) <- list(Gnames,c("lower","upper"), paste("E",1:6,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.95)[1:2]
  }
}

#Figure A2
par(mfrow=c(1,2))

for (k in 1:2){  
plot(1:m,MCMC.covtensor$av.G.coord[,k,],ylab="",xlab="",pch=16,xaxt="n",frame.plot=F,xlim=c(0.5,m+.5),ylim=c(-.3,.5),main = "")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,MCMC.covtensor$av.G.coord[,k,],1:m,HPD.tensor.coord[,1,k],length=0.1,angle=90)
arrows(1:m,MCMC.covtensor$av.G.coord[,k,],1:m,HPD.tensor.coord[,2,k],length=0.1,angle=90)
mtext(dimnames(MCMC.covtensor$av.G.coord)[[2]][k],side=3,at=0,font=2)
}

round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)


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

cor(e11,eigen(gamma)$vectors)
cor(e12,eigen(gamma)$vectors)
cor(e13,eigen(gamma)$vectors)

cor(e21,eigen(gamma)$vectors)
cor(e22,eigen(gamma)$vectors)
cor(e23,eigen(gamma)$vectors)

par(mfrow=c(2,2))
par(mar=c(5,4,4,2))
plot(1:m,rowMeans(e11.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e11))),las=2)
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,rowMeans(e11.proj),1:m,HPD.e11[,1],length=0.1,angle=90)
arrows(1:m,rowMeans(e11.proj),1:m,HPD.e11[,2],length=0.1,angle=90)
mtext("e11",side=3,at=0,font=2)

plot(1:m,rowMeans(e12.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e11))),las=2)
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,rowMeans(e12.proj),1:m,HPD.e12[,1],length=0.1,angle=90)
arrows(1:m,rowMeans(e12.proj),1:m,HPD.e12[,2],length=0.1,angle=90)
mtext("e12",side=3,at=0,font=2)

plot(1:m,rowMeans(e21.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e21))))
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,rowMeans(e21.proj),1:m,HPD.e21[,1],length=0.1,angle=90)
arrows(1:m,rowMeans(e21.proj),1:m,HPD.e21[,2],length=0.1,angle=90)
mtext("e21",side=3,at=0,font=2)

plot(1:m,rowMeans(e22.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e22))))
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,rowMeans(e22.proj),1:m,HPD.e22[,1],length=0.1,angle=90)
arrows(1:m,rowMeans(e22.proj),1:m,HPD.e22[,2],length=0.1,angle=90)
mtext("e22",side=3,at=0,font=2)

save(list=ls(),file='~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_tensor.RData')


