# utility functions for Gene-level quantitative trait mapping in Caenorhabditis elegans.
# Luke M. Noble, Matthew V. Rockman, Henrique Teotónio
# June 2019

simPheno <- function(h2, Ncausal, X, returnCausalIx=FALSE){
  
  # create phenoptypes with a specified heritability.
  # h2 is vector of additive and epistatic (pairwise interaction) heritabilities
  # from Normal effects.
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
  
  # sample additive effects
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

assignRecDoms <- function(df, recf='rec.domains') {
  #ordered df of 4 breakpoints per chromosome
  recd <- fread(recf, header=T, data.table = F)[,1:5]
  recd <- melt(recd, id='chrom', value.name='posr', variable.name = 'junction')
  recd$posr <- recd$posr*1e6
  stopifnot(nrow(recd)==24)
  df$domain = NA
  do.call('rbind', lapply(split(df, df$chrom), function(x) {
    rec = recd[recd$chrom==x$chrom[1],]
    x$domain[x$pos <= rec$posr[1]] <- 'tip'
    x$domain[is.na(x$domain) & x$pos <= rec$posr[2]] <- 'arm'
    x$domain[is.na(x$domain) & x$pos <= rec$posr[3]] <- 'center'
    x$domain[is.na(x$domain) & x$pos <= rec$posr[4]] <- 'arm'
    x$domain[is.na(x$domain)] <- 'tip'
    x
  }))
}

LDprune <- function(snps, x, maxr2 = 0.99, window=3000, step=1000){
  # x=marker x individual matrix
  P = nrow(x)
  X = t(x)
  wins = seq(1, P, window)
  wins = c(sort(c(wins, seq(step, P, window))), P+1)
  keepers = outs = NULL
  for(i in 1:(length(wins)-1)){
    ix = wins[i]:(wins[i+1]-1)
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

contiguousRL <- function(x, fillRun=NULL, flattenRun=NULL){
  # assign unique integer to contiguous runs of x
  # optionally force |distances| < flattenRun to 1
  # and fill runs of <= minRun, flanked by contiguous matching runs of > minRun
  if(!is.null(flattenRun)) x[abs(x)>1 & abs(x)<flattenRun] <- 1
  rl = rle(x)
  if(!is.null(fillRun)){
    i = which(rl$lengths <= fillRun)
    i = i[i>1 & i<len(rl$lengths)]
    ll = rl$lengths[i-1]
    lr = rl$lengths[i+1]
    ok = i[ll > fillRun & lr > fillRun & rl$values[i-1]==rl$values[i+1]]
    if(len(ok)>0){
      rl$lengths[ok-1] = rl$lengths[ok-1]+rl$lengths[ok]
      rl$lengths = rl$lengths[-ok]  
    }
  }
  rep.int(1:len(rl$lengths), times = rl$lengths)
}
