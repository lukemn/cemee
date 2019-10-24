#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)
Xfile1 = args[1]              # pruned marker X line genotype matrix for 1 chromosome (.h5 file, in slot 'X')
Xfile2 = args[2]              # pruned marker X line genotype matrix for 1 chromosome (.h5 file, in slot 'X')
ha = as.numeric(args[3])      # additive heritability
haa = as.numeric(args[4])     # pairwise epistatic heritability
Ncausal = as.integer(args[5]) # number of causal markers per simulation
Nsim = as.integer(args[6])    # number of simulations to run
prefix = args[7]              # output prefix (e.g., job array ID)
require(data.table, quietly = T, warn.conflicts = F)
require(rhdf5, quietly = T)

# write Nsim simulated phenotypes for heritability (ha, haa), Ncausal markers to file 
# 1. genotypes and trait values: sim_pref_i_ha_haa_Ncausal.phe.csv for i in 1:Nsim (Nlines x Ncausal + y)
# 2. sampled positions : sim_pref_ha_haa_Ncausal.ix.csv (Nsim x Ncasual)

X = rbind(h5read(Xfile1, 'X'), h5read(Xfile2, 'X'))

simPheno <- function(h2, Ncausal, X, returnCausalIx=FALSE){
  
  # generate phenotypes with a specified heritability, effect sizes from N(0,1).
  # h2 is vector of additive and epistatic (pairwise interaction) heritabilities
  # nCausal is number of causal variants to be sampled from X,
  # the M marker x N line genotype matrix. 
  # Each heritability component will be based on the same set of Ncausal markers.
  # output is data.frame of x1:xNcausal matrix and y=phenotype vector
  # or list of [[sampled positions]][[data.frame]] if returnCausalIx=TRUE.
  # In the special case of 2-marker epistasis, genotypes are filtered to ensure
  # 1. presence of all four genotype classes at a minimum frequency of 3
  # 2. a maximum of 5 NA/ambiguous [%%1 > 0] genotypes
  
  M <- dim(X)[1]
  N = dim(X)[2]
  stopifnot(all(M>1, N>1, M>Ncausal, length(h2)==2, sum(h2)<1, sum(h2)>0))
  minClassFreq = 3
  maxNA = 5
  
  # sample causal loci
  if(h2[1]>0 & Ncausal==2){
    while(1){
      ix_causal <- sample(M, Ncausal)
      x = M[ix_causal,]
      xtab = table(x)
      if(min(dim(xtab))==2 & min(xtab)>=minClassFreq & sum(unlist(x)%%1>0 | is.na(xtab))<=maxNA) break
    }
  } else {
    ix_causal <- sample(M, Ncausal)
  }
  # save for later
  Xraw <- X
  # scale sampled genotypes to mean 0, stdev 1
  X <- as.matrix(X[ix_causal,])
  X <- scale(t(X))
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
  if(Ncausal>1) x = t(x)
  df = data.table(x, y = Y)
  names(df)[1:Ncausal] = paste0('x', 1:Ncausal)
  if(returnCausalIx) return(list(ix_causal, df))
  return(df)
}

sims = lapply(1:Nsim, function(i) simPheno(c(ha, haa), Ncausal, X, returnCausalIx = T))
ix = data.table(do.call(rbind, lapply(sims, '[[', 1)))
sims = lapply(sims, '[[', 2)
h5h = sprintf('sim_%s_%s_%s_%s.h5', prefix, ha, haa, Ncausal)
h5createFile(h5h)
h5write(sims, h5h, 'phe')
h5write(ix, h5h, 'ix')
