

simPheno <- function(h2, Ncausal, X, returnCausalIx=FALSE){
  
  # friendlier version of simY, A & AA only
  # this function creates phenoptypes with a specified heritability.
  # h2 is vector of additive and epistatic (pairwise interaction) heritabilities
  # nCausal is number of causal variants to be sampled from X
  # X is the M marker x N line genotype matrix
  # output is data.frame of x1:xNcausal matrix and y=phenotype vector
  # or list of [[sampled positions]][[data.frame]] if returnCausalIx=TRUE
  
  M <- dim(X)[1]
  N = dim(X)[2]
  stopifnot(all(M>1, N>1, M>Ncausal, length(h2)==2, sum(h2)<1, sum(h2)>0))
  
  # sample causal loci
  ix_causal <- sample(M, Ncausal)
  Xraw <- X
  # scale genotypes to mean 0, stdev 1
  X <- scale(t(X))
  X <- as.matrix(X[,ix_causal])
  # mean impute missing data
  for(i in 1:Ncausal) X[,i][is.na(X[,i])] <- mean(X[,i], na.rm=T)
  Mc <- Ncausal
  h2 <- c(h2, 1-sum(h2))
  
  # sample additive effects from the standard normal distribution
  effA <- rnorm(Mc)
  YA <- X %*% effA
  YA <- YA %*% sqrt(h2[1]/var(YA))
  Y = YA
  
  # sample epistatic effects
  if(Mc>1 & h2[2]>0){
    M2d = (Mc*(Mc-1))/2
    XAA <- matrix(0, N, M2d)
    for (i in 1:N){
      gaa <- X[i,] %*% t(X[i,])
      XAA[i,] <- gaa[upper.tri(gaa)]
    }  
    effAA <- as.matrix(rnorm(M2d))
    Y2d <- XAA %*% effAA
    Y2d <- Y2d %*% sqrt(h2[2]/var(Y2d))
    Y = Y + Y2d
  }
  
  # add residual
  E <- rnorm(N)
  E <- E * sqrt(h2[3]/var(E))
  Y = as.numeric(scale(Y+E))
  
  x = Xraw[ix_causal,]
  df = cbind(t(x), data.frame(y = Y))
  names(df)[1:Ncausal] = paste0('x', 1:Ncausal)
  if(returnCausalIx) return(list(ix_causal, df))
  return(df)
}