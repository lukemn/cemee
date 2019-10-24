#!/usr/bin/env Rscript

# Estimation of LD effects on G estimated from additive genetic similarity,
# by varying within population pairwise marker pruning.

require(sommer)
require(parallel)
require(data.table)

load('ngm6_G_transitions.rda', verbose=T)         # plate level transition rate estimates corrected for technical covariates
gt <- fread('S1_WS220_CeMEEv2_markerSet1.csv.gz') # genotypes

NP=4     # parallel threads
nrep=10  # number of iterations of random pruning (within each window, a single marker from those in LD > r2 threshold is chosen randomly each iteration)
r2s = c(1, 0.99, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5) # LD pruning thresholds

# Lines per population. 
# For A6140, lines will be drawn from the CA population sizes each iteration.
popn = table(unique(ngm[,c('line', 'poprep')])$poprep)

# Extract G/E, additive genetic correlations and heritability, 
# breeding values, fixed effect coefficients, and residuals from a sommer model.
# random effect model including id/units
mmsum <- function(model){
  
  reff = names(model$sigma_scaled)
  g <- model$sigma_scaled[[grep('id', reff)]]
  e <- model$sigma_scaled[[grep('units', reff)]]
  
  list(G=g,
       E=e,
       gen.cor = cov2cor(g),
       h2 = diag(g) / (diag(g) + diag(e)),
       bv = model$U[[grep('id', reff)]],
       fixeff = model$Beta,
       resid = model$residuals
  )
}

# Additive genetic similarity matrix by SVD on standardised genotypes (gmat = P markers x N lines)
hgsm <- function(gmat){
  m = nrow(gmat) # markers
  gmat = scale(t(gmat))
  gmat[is.na(gmat)] <- 0
  # u = PCs, v = marker loadings, d = singular values for each PC
  # so that xx = svdx$u %*% diag(svdx$d) %*% t(svdx$v). K = XX'= US(US)'
  X = tcrossprod(gmat) / (m-1)
  svdx <- svd(X)
  D=diag(sqrt(svdx$d))
  US = svdx$u %*% D
  K = tcrossprod(US)
  K <- K/mean(diag(K))
  rownames(K) <- colnames(K) <- rownames(gmat)
  K
}

# Pairwise LD prune (x = P x N matrix for a single chromosome)
LDprune <- function(snps, x, maxr2 = 0.99, window=3000, step=2000, randomSeed=F){
  
  P = nrow(x)
  X = t(x)
  if(window > P) {
    wins = c(1, P+1)
  } else {
    wins = seq(1, P, window)
    wins = c(sort(c(wins, seq(step, P, window))), P+1)  
  }
  keepers = outs = NULL
  for(i in 1:(length(wins)-1)){
    ix = wins[i]:(wins[i+1]-1)
    if(randomSeed) ix = sample(ix)
    rx = cor(X[,ix])^2
    rx[upper.tri(rx, diag = T)] <- 0
    for(j in 1:nrow(rx)){
      if(!j %in% outs){
        jj = ix[which(rx[,j]>maxr2)]
        outs = append(outs, jj)
      }
    }
    outs = unique(outs)
    keepers = unique(append(keepers, ix[!ix %in% outs]))
  }
  list(snps[keepers,], x[keepers,])
}


mclapply(split(ngm, ngm$poprep), mc.cores=NP, function(x) {
  
  lapply(1:nrep, function(k) {

    pop = x$poprep[1]
    x$date = factor(x$date) # block 
    if(pop=='A6140') x = subset(x, line %in% sample(unique(x$line), sample(popn[-1], 1)))
    
    for(r2 in r2s){
      
      o=sprintf('mfit_%s_r2_%s_%s_sum.rda', pop, k, r2)
      if(!file.exists(o)){
        snps = as.data.frame(gt[,1:2])
        X = as.data.frame(gt[,colnames(gt) %in% x$line,with=F])
        X[X==1] <- 2
        X[X==0] <- 0
        X[X>0 & X<1] <- 1 # intermediate HMM genotypes are set to het
        # limit to segregating markers
        seg <- (apply(X, 1, sum)/(ncol(X)*2)) %% 1 > 0
        X <- X[seg,]
        snps = snps[seg,]
        # sequential LD prune
        if(r2<1){
          print(table(snps$chrom))          
          while(1){
            i = nrow(snps)
            lds <- mclapply(split(cbind(snps, X), snps$chrom), mc.cores = NP, function(i) LDprune(i[,1:2], i[,-(1:2)], window=as.integer(1500*r2), step=as.integer(1000*r2), maxr2 = r2, randomSeed=T))
            snps = do.call(rbind, lapply(lds,  '[[', 1))
            X = do.call(rbind, lapply(lds,  '[[', 2))
            if(nrow(snps)==i) {break; print(table(snps$chrom))}
          }
        }
        
        # Make the GSM, fit the model
        Ax <- hgsm(X)
        mfit = mmer(cbind(sf, sb, fs, fb, bs, bf) ~ 1, random = ~ vs(id, Gu=Ax), rcov=~units, data=x)
        msum <- mmsum(mfit)
        save(msum, file = o)
      }
    }  
  })
})

