rm(list = ls())
library(data.table)
library(pracma)
library(matrixStats)
library(ggplot2)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

# Load the WI matrices
N2_mat=as.matrix(read.table("~/PATH/TO/DIR/MA_lines/VM_estimates_N2.txt",h=TRUE,sep="\t"))/5
PB_mat= as.matrix(read.table("~/PATH/TO/DIR/MA_lines/VM_estimates_PB.txt",h=TRUE,sep="\t"))/5

lmat=cbind(0,1,2,3)
layout(lmat,widths=c(1,rep(5,3)), heights=c(5))

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}

#quartz()

#### Angle for the first 3 EV

i1=1
par(yaxt="s",mar=c(0.2,2,0,0.2),xaxt="n")
for(i1 in 1:3){
A = VCV_mat[[1]]$G1_mat/2
eigVal_A <- eigen(A)$values; eigVec_A <- eigen(A)$vectors
yvect_gmax=cbind()

plot(c(-1, 1), c(0, 0), type = "b", lwd = 2, ylim = c(-1,1), asp = 1, 
		ylab = "", yaxt = "n",bty="n",las=1,xlab="")

axis(side=2,at=c(-1,0,1),labels=c("","",""),pos=-1.1)

if(i1==1){
	mtext(side=2,expression(atop(paste("Orientation (",theta,") and"))),cex=1.5,padj=-.2)
	mtext(side=2,expression(paste("\n magnitude of ",g[i])),padj=-.6,cex=1.5)
}

for(temp_c in 2:7){

temp_G <- VCV_mat[[temp_c]]$G1_mat/2
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])
temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))
		points(c(-temp_coordinates[1], temp_coordinates[1]), c(-temp_coordinates[2], temp_coordinates[2]), type = "l", col = rep(c("dodgerblue","firebrick3"),c(4,3))[temp_c], lwd = 2)
}


temp_G <- N2_mat
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])
temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))
points(c(-temp_coordinates[1], temp_coordinates[1]), c(-temp_coordinates[2], temp_coordinates[2]), type = "l", col ="darkgreen", lwd = 3,lty=2)

temp_G <- PB_mat
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])
temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))
points(c(-temp_coordinates[1], temp_coordinates[1]), c(-temp_coordinates[2], temp_coordinates[2]), type = "l", col ="magenta", lwd = 3,lty=2)

if(i1==1) mtext(side=1,expression(paste(g[1]," (62%) ")),padj=-.7,cex=1.5)
if(i1==2) mtext(side=1,expression(paste(g[2]," (23%) ")),padj=-.7,cex=1.5)
if(i1==3) mtext(side=1,expression(paste(g[3]," (8%) ")),padj=-.7,cex=1.5)

}


