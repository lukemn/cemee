rm(list = ls())
library(data.table)
library(pracma)
library(matrixStats)
library(ggplot2)


# Produce Fig 4


angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}
load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData')

# Load the simulation results
Simul = new.env()
load("~/PATH/TO/DIR/Simulations_AmNat/Main_model_output_list_with_metadata.RData", envir = Simul)

lmat=t(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13),4,4))
layout(mat=lmat, heights=c(1,1,1,.7),widths=c(.4,1,1,1))
### first plot simply for the Y axis
par(mar=c(5,4,4,2))
plot(0,0,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")

mtext(side=2,expression(atop(paste("Orientation (",theta,") and magnitude of ",g[i]))),adj=.7,padj=1)

par(mar=c(0,0,0,0))
x_vect=seq(0,1,.01)
y_vect=sqrt(1-x_vect^2)


for(kkk in 1:3){
	if(kkk==1) vect_j = 2:4
	if(kkk==2) vect_j = 5:7	
	if(kkk==3) vect_j = 1:20

#### Angle for the first 3 EV
#i1=1
#par(yaxt="s",mar=c(0.2,0.7,0,0.2),xaxt="n")
par(yaxt="s",xaxt="n")
for(i1 in 1:3){
A = VCV_mat[[1]]$G1_mat/2
eigVal_A <- eigen(A)$values; eigVec_A <- eigen(A)$vectors
yvect_gmax=cbind()

plot(c(0, 1), c(0, 0), type = "b", lwd = 2, ylim = c(-.2,1.2), asp = 1, 
		ylab = "", yaxt = "n",bty="n",las=1,xlab="")

if(kkk==1 & i1==1) text(0.1,1.15,"CA[1-3]50",cex=2)
if(kkk==2 & i1==1) text(0.2,1.13,"CA[1-3]100",cex=2)
if(kkk==3 & i1==1) text(0.2,1.15,"Simulated G",cex=2) 

axis(side=2,at=c(0,1),labels=c("",""),pos=-0.1)

### We will subsample and create an area for each CA[1-3]G50

if(kkk!=3){
vect_sampling=sample(1:nrow(VCV_mat[[1]]$VCV_Mat),1000)

list_sampled_coordinates=list()
k_count=0
for(j in vect_j){
	temp_coordinates=NULL
for(ks in vect_sampling){
	
	A_temp = matrix(VCV_mat[[1]]$VCV_Mat[ks,1:36],6,6)/2
	temp_G <- matrix(VCV_mat[[j]]$VCV_Mat[ks,1:36],6,6)/2

	temp_angle <- angle_eigenV(eigen(A_temp)$vectors[, i1], eigen(temp_G)$vectors[, i1])

temp_coordinates <- rbind(temp_coordinates ,c(temp_angle, eigen(temp_G)$values[i1]/eigen(A_temp)$values[i1]))
}
k_count=k_count+1
list_sampled_coordinates[[k_count]]= temp_coordinates
}
}
k=0
for(j in vect_j){

	k=k+1

if(kkk!=3){
temp_G <- VCV_mat[[j]]$G1_mat/2
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])
if(temp_angle>(pi/2)) temp_angle= pi-temp_angle
temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))

temp_sampled_angles=list_sampled_coordinates[[k]][,1]
temp_sampled_int=list_sampled_coordinates[[k]][,2]

temp_sampled_angles[temp_sampled_angles > (pi/2)]=pi-temp_sampled_angles[temp_sampled_angles > (pi/2)]

angles_for_arrow=sort(temp_sampled_angles)[c(0.05,0.95)*length(vect_sampling)]
int_for_arrow=sort(temp_sampled_int)[c(0.05,0.95)*length(vect_sampling)]

new_coord_for_arrows=rbind(pol2cart(c(angles_for_arrow[1], eigen(temp_G)$values[i1]/eigVal_A[i1])),pol2cart(c(angles_for_arrow[2], eigen(temp_G)$values[i1]/eigVal_A[i1])))

vect_col =c("black","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick","forestgreen")

arrows(new_coord_for_arrows[1,1],new_coord_for_arrows[1,2],new_coord_for_arrows[2,1],new_coord_for_arrows[2,2],code=3,length=.1,angle=90,col=vect_col[j])

new_coord_for_arrows=rbind(pol2cart(c(temp_angle, int_for_arrow[1])),pol2cart(c(temp_angle, int_for_arrow[2])))

arrows(new_coord_for_arrows[1,1],new_coord_for_arrows[1,2],new_coord_for_arrows[2,1],new_coord_for_arrows[2,2],code=3,length=.1,angle=90,col=vect_col[j])
points(c(0, temp_coordinates[1]), c(0, temp_coordinates[2]), type = "l", col = vect_col[j], lwd = 4)

}else{
	vect_col = c("green")
	temp_G <- Simul$G_mat_final[[j]]/2
temp_angle <- angle_eigenV(eigen(A)$vectors[, i1], eigen(temp_G)$vectors[, i1])
if(temp_angle>(pi/2)) temp_angle= pi-temp_angle
temp_coordinates <- pol2cart(c(temp_angle, eigen(temp_G)$values[i1]/eigVal_A[i1]))
points(c(0, temp_coordinates[1]), c(0, temp_coordinates[2]), type = "l", col = vect_col, lwd = 2)

	}

}

lines(x_vect, y_vect,lty=2)
}}

par(mar=c(5,4,4,2))
plot(0,0,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
mtext(side=1,expression(paste(g[1]," (62%) ")),padj=-5,cex.axis=1.8)
plot(0,0,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
mtext(side=1,expression(paste(g[2]," (23%) ")),padj=-5,cex.axis=1.8)
plot(0,0,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
mtext(side=1,expression(paste(g[3]," (8%) ")),padj=-5,cex.axis=1.8)








