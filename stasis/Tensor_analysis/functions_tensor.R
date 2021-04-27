
inv.rootP <- function(P) {
	rootP <- matrix(0, n, n)
	for (i in 1:n) {
		val <- eigen(P)$values
		vec <- eigen(P)$vectors
		rootP <- rootP + (vec[, i] %*% t(vec[, i])) * sqrt(val[i])
	}
	solve(rootP)
}


##### covtensor function
#START
library(matrixcalc)
library(gdata)

covtensor <- function(Gs) {
	if (dim(Gs)[[1]] != dim(Gs)[[2]]) {
		stop("G array must be of order n x n x m x MCMCsamp")
	}
	if (is.na(dim(Gs)[4])) {
		stop("There are no MCMCsamples")
	}
	neigten <- n * (n + 1)/2
	#Number of eigentensors
	MCMC.S <- array(, c(neigten, neigten, MCMCsamp))
	dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep = ""), paste("e", 1:neigten, sep = ""))
	for (k in 1:MCMCsamp) {
		MCMCG <- Gs[, , , k]
		MCMCvarmat <- t(apply(MCMCG, 3, diag))
		#find the variances of the kth G and store them
		MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle))
		#find the covariances of the kth G and store them
		MCMC.S[1:n, 1:n, k] <- cov(MCMCvarmat, MCMCvarmat)
		#fill the upper left quadrant of the kth S
		MCMC.S[(n + 1):neigten, (n + 1):neigten, k] <- 2 * cov(MCMCcovmat, MCMCcovmat)
		#fill the lower right quadrant of the kth S
		MCMC.S[1:n, (n + 1):neigten, k] <- sqrt(2) * cov(MCMCvarmat, MCMCcovmat)
		#fill the upper right quadrant of the kth S
		MCMC.S[(n + 1):neigten, 1:n, k] <- sqrt(2) * cov(MCMCcovmat, MCMCvarmat)
		#fill the lower left quadrant of the kthS
		}
	av.S <- apply(MCMC.S, 1:2, mean)
	#posterior mean S
	av.S.val <- eigen(av.S)$values
	#eigenvalues of posterior mean S
	av.S.vec <- eigen(av.S)$vectors
	#eigenvalues of posterior mean S
	eTmat <- array(, c(n, n, neigten))
	dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep = ""))
	for (i in 1:neigten) {
		emat <- matrix(0, n, n)
		lowerTriangle(emat) <- 1/sqrt(2) * av.S.vec[(n + 1):neigten, i]
		emat <- emat + t(emat)
		diag(emat) <- av.S.vec[1:n, i]
		eTmat[, , i] <- emat
	}
	#construct the second-order eigentensors of posterior mean S
	eT.eigen <- array(, c(n + 1, n, neigten))
	for (i in 1:neigten) {
		eT.eigen[1, , i] <- t(eigen(eTmat[, , i])$values)
		#Eigenvalues of the ith eigentensor
		eT.eigen[2:(n + 1), , i] <- eigen(eTmat[, , i])$vectors
		#Eigenvectors of the ith eigentensor
		eT.eigen[, , i] <- eT.eigen[, order(abs(eT.eigen[1, , i]), decreasing = T), i]
	}
	MCMC.S.val <- matrix(, MCMCsamp, neigten)
	colnames(MCMC.S.val) <- paste("E", 1:neigten, sep = "")
	for (i in 1:MCMCsamp) {
		for (j in 1:neigten) {
			MCMC.S.val[i, j] <- t(av.S.vec[, j]) %*% MCMC.S[, , i] %*% av.S.vec[, j]
		}
	}
	#posterior distribution of the genetic variance for the eigenvectors of posterior mean S
	av.G.coord <- array(, c(m, neigten, 1))
	dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep = ""))
	for (i in 1:neigten) {
		av.G.coord[, i, ] <- apply((apply(Gs, 1:3, mean)), 3, frobenius.prod, y = eTmat[, , 
			i])
	}
	#Coordinates of the jth avG for the eigentensors of posterior mean S
	MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
	dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep = ""))
	for (i in 1:neigten) {
		MCMC.G.coord[, i, ] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[, , i])
	}
	#Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
	tensor.summary <- data.frame(rep(av.S.val, each = n), t(data.frame(eT.eigen)))
	colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
	rownames(tensor.summary) <- paste(paste("e", rep(1:neigten, each = n), sep = ""), rep(1:n, 
		neigten), sep = ".")
	list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, 
		MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}
#END

# Rotate matrix 270 clockworks
flip.matrix <- function(x) {
	mirror.matrix(rotate180.matrix(x))
}
# Mirror matrix (left-right)
mirror.matrix <- function(x) {
	xx <- as.data.frame(x)
	xx <- rev(xx)
	xx <- as.matrix(xx)
	xx
}
# Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
	mirror.matrix(t(x))
}


R.proj <- function(Gs,p,vec){
 if (dim(Gs)[[1]] != dim(Gs)[[2]]){
   stop("G array must be of order n x n x m x MCMCsamp")
 }
 if (is.na(dim(Gs)[4])) {
   stop("There are no MCMCsamples")
 }
 n <- dim(Gs)[[1]]
 m <- dim(Gs)[[3]]
 MCMCsamp <- dim(Gs)[[4]]
  rand.vec <-matrix(,vec,n)
    for (i in 1:vec){
      b <- runif(n,-1,1)
      rand.vec[i,] <- b/(sqrt(sum(b^2)))
    }
    #generate unit length random vectors  
    proj<- function(G,b) t(b) %*% G %*% (b)
    #internal function to do projection
      G.proj <- array(,c(MCMCsamp, m, vec))
        colnames(G.proj) <- dimnames(Gs)[[3]]
        for (i in 1:vec){
          G.proj[,,i]<- t(apply(Gs, 3:4, proj, b = rand.vec[i,]))
        }
      #project each random vector through each MCMC sample of each G
        prs <- cbind(rep(1:m, each = m), 1:m) 
        prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE] 
        #setting up an index for HPD comparisons
        proj.score <-matrix(,vec,((m^2 - m)/2))
          for (k in 1:vec){
            HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
            proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0) 
          }
          #for a given random vector, examine if the HPD intervals of any pair of G matrices overlap
      vec.score <-cbind(rand.vec, proj.score)
        colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
      #collate the random vectors and the outcome of their projection on the G matrices
      sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0) 
      #collate just the random vectors that resulted in significant differences in variance
        if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
        else{
          eig.R <- eigen(cov(sig.vec[,1:n]))
            rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
            colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
        }  
    #eigen analysis of the R matrix
    list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
}


proj<- function(G, b) t(b) %*% G %*% (b)



