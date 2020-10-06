#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)
prefix = as.integer(args[1])  # output prefix (e.g., job array ID, simulate if >0, summarise and plot if==0)
np = as.integer(args[2])      # parallel threads

require(data.table, quietly = T, warn.conflicts = F)
require(R.utils, quietly = T, warn.conflicts = F)
require(GridLMM, quietly = T)
require(parallel, quietly = T)
require(ggplot2, quietly = T)
require(Hmisc, quietly = T)
source('CeMEE_util.R')

############################################################
# focal SNP h2
h2fgs = 1:15/100
# number of SNPs in background component
nbgs = c(50, 100, 500)
# significance threshold (Multitrans)
thresh=3.23e-6
# window size for QTL interval (merged across qwindow/2 step)
qwindow=3e5
# QTL confidence interval
LODdrop=1.5
############################################################

## prepare the input files first, before parallel simulation.
prepD <- function(ofile = 'lmm_sim.rda'){
  
  ## Required in working directory:
  # WS220_CeMEEv2_markerSet1.csv.gz
  # WS220_CeMEEv2_markerSet2.csv.gz 
  # WS220_CeMEEv2_geneticMap.csv.gz
  # K.txt.gz (NxN genetic similarity matrix)
  
  # filter genotypes (<=2% hets, >5% MAF, r2 < 0.99)
  gtf <- fread('WS220_CeMEEv2_markerSet1.csv.gz')
  gtk <- fread('WS220_CeMEEv2_markerSet2.csv.gz')
  gmap <- fread('WS220.clean.geneticMap.csv.gz')
  hetkc <- apply(gtk[,-(1:4)], 1, function(x) sum(x %% 1 != 0))
  gtff <- gtf[hetkc<=(733*.02),]; dim(gtff)
  gtff <- filterMAF(gtff, MAFgt=0.05)
  snpi = gtff[,1:2]
  Xi = gtff[,-(1:2)]
  
  while(1){
    nin = nrow(snpi)
    print(sprintf('iter%s', i))  
    print(table(snpi$chrom))  
    lds <- parallel::mclapply(split(cbind(snpi, Xi), snpi$chrom), mc.cores = 6, function(i) LDprune(i[,1:2], i[,-(1:2)], window=2000, step=1000, maxr2 = 0.99))
    snpi = do.call(rbind, lapply(lds,  '[[', 1))
    Xi = do.call(rbind, lapply(lds,  '[[', 2))
    if(nrow(snpi)==nin){
      rownames(snpi) = rownames(Xi) = NULL  
      print(table(snpi$chrom))  
      break
    }
  }
  parallel::mclapply(split(Xi, snpi$chrom), mc.cores = 6, function(x) {x = as.matrix(x); summary(unlist(lapply(1:(nrow(x)-1), function(i) cor(x[i,], x[i+1,])^2)))})
  
  Xo = t(Xi)
  snpo = snpi
  rownames(snpo) = NULL
  colnames(Xo) = snpo$snp = paste0(snpo$chrom, ':', snpo$pos)
  K = as.matrix(fread('K.txt.gz'))
  rownames(K) = rownames(X)
  save(Xo, snpo, K, file = ofile)
  
  # write some sample phenotypes for multitrans using full foreground grid at nbg=100, h2=0.5
  sims = lapply(h2fgs, function(h2) simPolygenicPheno(c(h2, 0), c(0.5-h2, 0), 1, 100, Xt, returnCausalIx = F))
  o = data.frame(line=rownames(sims[[1]]), do.call(cbind, lapply(sims, '[', , 2)))
  fwrite(o, file = 'sim.phe', row.names = F, sep = '\t')
}

## main simulation function
simulateQ <- function(fixedh2 = 0.5, Nsim = 200){
  
  load('lmm_sim.rda', verbose=F)
  pgrid = expand.grid(h2fgs, nbgs)
  names(pgrid) = c('h2', 'nbg')
  pgrid$i = rownames(pgrid)
  
  apply(pgrid, 1, function(i) {
    ii = as.numeric(i)
    h2 = ii[1]; nbg = ii[2]
    h2bg = fixedh2-h2
    out = do.call(rbind, mclapply(1:Nsim, mc.cores=np, function(k) {
      simi = simPolygenicPheno(c(h2, 0), c(h2bg, 0), 1, nbg, Xt, returnCausalIx = T)
      df = simi[[1]]; fg = simi[[2]]; bg = simi[[3]]
      df$line = rownames(df)
      gwas = GridLMM_GWAS(formula = y~1+(1|line), test_formula = ~1, reduced_formula = ~0, 
                          data = df, X = Xo, X_ID = 'line', relmat = list(line = K), 
                          algorithm = 'Full', mc.cores = 1)
      res = cbind(snpo, p = round(-log10(gwas$results$p_value_REML), 3))
      res$causal = 0
      res$causal[rownames(res) %in% bg] = 1
      res$causal[rownames(res) %in% fg] = 2
      res$r2 = round(summary(lm(y~x1,df))$adj, 3)
      res
    }))
    out$h2 = h2
    out$nbg = nbg
    of = sprintf('sim_%s_%s.txt', prefix, i[3])
    fwrite(out, file = of, row.names = F, sep='\t')
    gzip(of)
  })
  
}

## define QTL and generate summary stats
crunchQ <- function(NP=20){
  
  # batch process simulations (each file is for a single scenario)
  
  simstat <- function(f, qwindow, LODdrop, thresh, minMarkerPerW=3, np=1){
   
    # get number and size of all detected QTL intervals
    # for any detected SNP, look within qwindow bp for 
    # other significant snps ranked by LD
    
    print(f)
    simb <- fread(f)
    simb$snp = paste0(simb$chrom, ':', simb$pos)
    nsim = nrow(simb)/nrow(snpo)
    ix = rep(1:nsim, each = nrow(snpo))
    
    mergeWindows <- function(wins, qwindow, minMarkerPerW){
      # merge windows with <minMarkerPerW markers, shrinking toward the center
      wt = table(wins)
      m = median(unique(wins))
      while(min(wt)<minMarkerPerW){
        w = wins[wins %in% names(wt)[wt<minMarkerPerW]]
        wmerge = w+qwindow*sign(m-w)
        wins[wins %in% names(wt)[wt<minMarkerPerW]] <- wmerge
        wt = table(wins)
      }
      wins
    }
    
    windowsToIntervals <- function(y, wins, LODdrop){
      
      # y is a data.frame of p-values for a single chromosome.
      # Given a vector of windows to split by (length = N snps),
      # return a bool vector indicating the QTL interval from 
      # the (single) peak marker (-log10(p) - LODdrop), calculated
      # from markers ranked by ld with the peak marker. 
      # Minimum interval is peak+-flanking.
      
      as.numeric(unlist(lapply(split(y, wins), function(z) {
        if(sum(z$det)==0) {
          rep(F, nrow(z))
        } else {
          pkz = which.max(z$p)
          # rerank by ld with peak marker
          zld = crossprod(Xs[,z$snp])
          z$ld = zld[,pkz]
          zo = z[order(sign(z$pos-z$pos[pkz]), zld[,pkz]),]
          pkz = which.max(zo$p)
          if(pkz<nrow(z) & pkz>1) zo[(pkz+1):nrow(z),] <- zo[nrow(z):(pkz+1),]
          pk = which.max(zo$p)
          qint = zo$p > (v$p[pk]-LODdrop)
          # set flanking (unless terminal, then nearest)
          if(sum(qint)<3) {
            l = ifelse(pk==1, 0, 1)
            r = ifelse(pk==nrow(zo), 0, 1)
            qint[(pk-l):(pk+r)] = T
          }
          qint = range(zo$pos[qint])
          o = rep(F, nrow(z))
          o[z$pos >= qint[1] & z$pos <= qint[2]] = T
          o
        }
      })))
    }
    
    # loop over simulations, return merged results
    o = do.call(rbind, mclapply(split(simb, ix), mc.cores = np, function(x){
      x$det = x$p > -log10(thresh)
      x$maf = mafs
      x$qtl = 0
      # loop over chromosomes
      x = do.call(rbind, lapply(split(x, x$chrom), function(y){
        if(sum(y$det)>0){
          # merge QTL intervals from stepped windows
          wins = mergeWindows(floor(y$pos/qwindow)*qwindow, qwindow, minMarkerPerW)
          q = windowsToIntervals(y, wins, LODdrop)
          wins = mergeWindows(ceiling((y$pos+qwindow/2)/qwindow)*qwindow, qwindow, minMarkerPerW)
          y$qtl = q | windowsToIntervals(y, wins, LODdrop)
          # QTL intervals cannot exceed termini of the merged window
          y$qtl[1] = 0
          y$qtl[length(y$qtl)] = 0
          y
        } else {
          y
        }
      })) # end chromosome
      x$qtlix = contiguousRL(x$qtl)
      x$qtlix[!x$qtl] = NA
      qw = table(x$qtlix)
      # size of focal QTL (if detected)
      simint = simq = pkd = x$qtlix[x$causal==2]
      if(!is.na(simq)) {
        i = subset(x, qtlix==simq)
        simint = diff(range(i$pos))
        # distance of simulated site from QTL peak (if detected)
        # take min of focal/background in rare cases where both fall within the focal 
        # QTL interval and the latter has a smaller pvalue
        pkd = abs(diff(c(i$pos[i$causal==2], i$pos[which.max(i$p)])))
        if(sum(i$causal==1)>0){
          if(max(i$p[i$causal==1]) > i$p[i$causal==2]){
            bgpk = which.min(abs(i$pos[i$causal==1] - i$pos[which.max(i$p)]))
            pkd = min(c(pkd, abs(diff(c(i$pos[bgpk], i$pos[which.max(i$p)])))))
          }
        }
      }
      # total number of qtl (exc. focal)
      nqtl = length(qw[names(qw)!=simq])
      nfp = nfp_exc = 0
      if(nqtl>0){
        # number of false positives (excluding focal and background QTL)
        qw = as.data.frame(table(x$causal, x$qtlix)>0)
        nfp = sum(qw[1,] & !apply(qw[-1,,drop=F], 2, any))  
        # on chromosomes other than the one with the focal QTL
        simchr = x$chrom[x$causal==2]
        qw = with(subset(x, chrom!=simchr), as.data.frame(table(causal, qtlix)>0))
        if(ncol(qw)>0) nfp_exc = sum(qw[1,] & !apply(qw[-1,,drop=F], 2, any))
      }
      cbind(x[x$causal==2,c('chrom', 'pos','p', 'det', 'r2', 'h2', 'nbg', 'maf')], interval=simint, nqtl, nfp, nfp_exc, pkd)
    })) # end simulation
    o
  }
  
  simf <- Sys.glob('sim*txt.gz')
  load('lmm_sim.rda', verbose=F)
  Xs = scale(Xo)
  mafs = round(afsToMafs(apply(Xo, 2, sum)/nrow(Xo)), 3)
  o = do.call(rbind, lapply(simf, function(f) simstat(f, thresh=thresh, qwindow=qwindow, LODdrop=LODdrop, np=NP)))
  save(o, file = sprintf('stat_w%s_lod%s_multrans.rda', qwindow, LODdrop))
}
  
# make some plots
plotQ <- function(qwindow, LODdrop){
  
  load(sprintf('stat_w%s_lod%s_multrans.rda', qwindow, LODdrop))
  o <- assignRecDoms(o)
  names(o)[names(o)=='domain'] <- 'Domain'
  
  # replace simulated h2 with true (estimated) h2
  o$sim_h2 = o$h2
  o$h2 = round(o$r2, 2)
  o = subset(o, h2 > 0 & h2 <= 0.2)
  o$mafq <- Hmisc::cut2(o$maf, g=4)
  
  # plotting range
  xl = c(0, 0.12)
  hsup = expression(italic(h^2))
  
  # where detected, but interval is NA, simulated QTL was outside the peak marker interval
  # check if the LODdrop interval is appropriately controlled
  qci = aggregate(data = subset(o, det), is.na(interval)~h2+nbg, function(x) sum(x)/len(x))
  names(qci)[3] = 'ci'
  qci$n = aggregate(data = subset(o, det), is.na(interval)~h2+nbg, len)[,3]
  # true site is within the interval ~99% of the time
  p <- ggplot(subset(qci, n>10), aes(h2, 1-ci)) + stat_summary() + theme_classic() + coord_cartesian(xlim = xl) + labs(x=hsup, y = 'Empirical confidence interval') + theme(plot.margin = margin(2, 10, 2, 2))
  ggsave('CI.png', h = 4, w = 4)
  
  # power by recombination rate domain
  p <- ggplot(o, aes(h2, as.numeric(det), col = Domain)) + stat_summary() + geom_smooth(method='glm', method.args = list(family='binomial'), se=F) + theme_classic() + coord_cartesian(xlim = xl) + labs(x=hsup, y = 'Power') + theme(legend.position = 'top', plot.margin = margin(2, 10, 2, 2))
  ggsave('power_by_domain.png', h = 4, w = 4)
  
  h2pred = data.frame(h2=rep(seq(2, 5, .01)/100, each=3))
  h2pred$Domain = rep(unique(o$Domain), times=(nrow(h2pred)/3))
  h2pred$power = predict(glm(det~as.numeric(h2)+Domain, o, family = 'binomial'), newdata = h2pred, type='resp')
  lint=.8
  melt(lapply(split(h2pred, h2pred$Domain), function(x) {x$h2[which.min(abs(x$power-lint))[1]]}))
  
  # pvals are higher for lower maf (for==r2, effect sizes are larger)
  # but the threshold is constant (inappropriately)
  p <- ggplot(subset(o, Domain!='tip'), aes(h2, as.numeric(det), col = mafq)) + stat_summary() + geom_smooth(method='glm', method.args = list(family='binomial'), se=F) + theme_classic() + coord_cartesian(xlim = xl) + labs(x=hsup, y = 'Power') + theme(legend.position = c(0.8, 0.3), plot.margin = margin(2, 10, 2, 2)) + scale_color_discrete('MAF quartile')
  ggsave('power_by_maf.png', h = 4, w = 4)
  
  # QTL interval by domain. centers are still large: median ~50-100kb
  aggregate(data = subset(o, det & h2 >= 0.03 & h2 <= 0.1), interval~h2+Domain, median)
  p <- ggplot(subset(o, det & Domain!='tip' & h2 >= 0.02), aes(h2, interval/1e3, col = Domain)) + stat_summary(fun.data = iqbox) + theme_classic() + coord_cartesian(xlim = xl) + labs(x=hsup, y = 'QTL interval (Kb)') + scale_y_log10() + theme(legend.position = 'top', plot.margin = margin(2, 10, 2, 2))
  ggsave('interval_by_domain.png', h = 4, w = 4)
  
  # distance of simulated focal marker from QTL peak marker
  aggregate(data = subset(o, det & h2 >= 0.03 & h2 <= 0.1), pkd~h2+Domain, median)
  aggregate(data = subset(o, det & h2 >= 0.03 & h2 <= 0.1), pkd~h2+Domain, function(x) quantile(x, .8))
  p <- ggplot(subset(o, det & Domain!='tip' & h2 >= 0.02), aes(h2, (pkd+1), col = Domain)) + stat_summary(fun.data = iqbox, position=position_dodge(width=0.005)) + theme_classic() + coord_cartesian(xlim = xl) + labs(x=hsup, y = 'Distance to peak marker (base pairs +1)') + scale_y_log10()  + theme(legend.position = 'top', plot.margin = margin(2, 10, 2, 2))
  ggsave('peak_distance_domain.png', h = 4, w = 4)
  
  # FDR (excluding chromosome with focal QTL)
  o$fdr = o$nfp_exc / o$nqtl
  o$fdr[o$det] = o$nfp_exc[o$det] / (o$nqtl[o$det]+1)
  p <- ggplot(subset(o, Domain!='tip' & h2>0), aes(h2, fdr, col = factor(nbg))) + stat_summary() + theme_classic() + coord_cartesian(xlim = xl) + labs(x=hsup, y = 'False discovery rate')  + theme(legend.position = 'top', plot.margin = margin(2, 10, 2, 2), legend.direction = 'vertical') + guides(color=guide_legend('Polygenicity', nrow=1))
  ggsave('fdr_nbg.png', h = 4, w = 4)
  
  # cf 2D power
  load('2dsim.ep3_mafq.rda', verbose=T)
  sub2d <- cbind(res[grep('tip', res$d, invert=T), c('h2', 'det', 'd')], test='2D')
  names(sub2d)[3] = 'Domain'
  pmerge <- rbind(sub2d, cbind(subset(o, Domain != 'tip')[,c('h2', 'det', 'Domain')], test='1D'))
  pmerge$test = factor(pmerge$test, levels = c('1D', '2D'))
  
  p <- ggplot(pmerge, aes(as.numeric(h2), as.numeric(as.logical(det)), group = Domain, col = Domain, linetype=test)) + stat_summary() + 
    theme_classic() + labs(y = 'Power', x=hsup) + theme(legend.position = 'right') + coord_cartesian(xlim = c(0, 0.1), ylim =c(0, 1.02)) + 
    guides(color = guide_legend('Domain(s)')) + scale_y_continuous(expand=c(0,0), breaks = seq(0, 1, .2)) + scale_x_continuous(expand=c(0,0)) +
    scale_shape('Model') + geom_smooth(method='glm', method.args = list(family = 'binomial'), se=F) + scale_linetype('Model')
  ggsave('power_1D2D_byDomain_nodrop.png', h=4, w=4.5)
}

# don't clobber!
if(!file.exists('lmm_sim.rda')) prepD()
if(prefix) simulateQ()
if(!prefix) {crunchQ(); plotQ()}
