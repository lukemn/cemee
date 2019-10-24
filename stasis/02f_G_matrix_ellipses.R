rm(list = ls())
library(data.table)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")

### Two plots that compose the Figure 3

quartz()

gen_CA = 100 # or set to 100
if(gen_CA==50){
vect_CA=2:4;ell_col="blue"	
}else{
vect_CA=5:7;ell_col="red"		
	}

lmat=rbind(
c(0,1:5),
c(6,11:15),
c(7,0,16:19),
c(8,0,0,20:22),
c(9,0,0,0,23:24),
c(10,0,0,0,0,25))
layout(lmat,w=c(.2,rep(1,5)),h=c(.2,rep(1,5)))

par(mar=c(0,0,0,0))
par( xaxt = "n", yaxt = "n", mar = c(.1, .1, .1, .1), bty = "n")
for(i in 2:6){ plot(1,1,type="n");text(1,1,c("SF","SB","FS","FB","BS","BF")[i],cex=1.5)}
for(i in 1:5){ plot(1,1,type="n");text(1,1,c("SF","SB","FS","FB","BS","BF")[i],cex=1.5)}



for (i1 in c(1:5)) {
	for (i2 in 2:6) {
		if (i1 < i2) {
			A = VCV_mat[[1]]$G1_mat[c(i1, i2), c(i1, i2)]
			ctr <- c(0, 0)
			angles <- seq(0, 2 * pi, length.out = 200)

			eigVal_A <- eigen(A)$values
			eigVec_A <- eigen(A)$vectors
			eigScl_A <- eigVec_A %*% diag(sqrt(eigVal_A)) # scale eigenvectors to length = square-root
			xMat_A <- rbind(ctr[1] + eigScl_A[1, ], ctr[1] - eigScl_A[1, ])
			yMat_A <- rbind(ctr[2] + eigScl_A[2, ], ctr[2] - eigScl_A[2, ])
			ellBase_A <- cbind(sqrt(eigVal_A[1]) * cos(angles), sqrt(eigVal_A[2]) * sin(angles)) # normal ellipse
			ellRot_A <- eigVec_A %*% t(ellBase_A) # rotated ellipse
			plot((ellRot_A + ctr)[1, ] * 1, (ellRot_A + ctr)[2, ] * 1, asp = 1, type = "n", 
				lwd = 2)

			lines((ellRot_A + ctr)[1, ], (ellRot_A + ctr)[2, ], asp = 1, type = "l", lwd = 1, lty = 1, col = "black")


			for (i in vect_CA) {
				B = VCV_mat[[i]]$G1_mat[c(i1, i2), c(i1, i2)]
				eigVal_B <- eigen(B)$values
				eigVec_B <- eigen(B)$vectors
				eigScl_B <- eigVec_B %*% diag(sqrt(eigVal_B)) # scale eigenvectors to length = square-root
				xMat_B <- rbind(ctr[1] + eigScl_B[1, ], ctr[1] - eigScl_B[1, ])
				yMat_B <- rbind(ctr[2] + eigScl_B[2, ], ctr[2] - eigScl_B[2, ])
				ellBase_B <- cbind(sqrt(eigVal_B[1]) * cos(angles), sqrt(eigVal_B[2]) * sin(angles)) # normal ellipse
				ellRot_B <- eigVec_B %*% t(ellBase_B) # rotated ellipse

				lines((ellRot_B + ctr)[1, ], (ellRot_B + ctr)[2, ], asp = 1, type = "l", lwd = 1, lty = 1, col = ell_col)

			}

	}
}
}

