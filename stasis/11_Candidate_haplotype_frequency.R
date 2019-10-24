#!/usr/bin/env Rscript

require(dplyr)

# haplotype statistics around candidate N2 lab adaptation alleles

# load gene positions, haplotype reconstructions filtered to 2Mb around each
load('~/Documents/cemee/elife/Github_Rcodes_clean/11_Candidate_haplotype_frequency_data.rda', verbose=T)
print(goi)
haplos = haploi

# mean haplotype length (cM) split by wild/N2
haplos$line = as.character(haplos$line)
haplos$haplotype = as.character(haplos$haplotype)
poph = substr(haplos$line, 1, 5)
mhlw <- do.call(rbind, lapply(split(haplos, poph), function(x) {
  adjUnit = 100
  cbind(pop = gsub('L', '', substr(x$line[1], 1, 6)),
        do.call(rbind, lapply(split(x, x$chrom), function(y){
          y$haplotype[y$haplotype %in% c('N2anc', 'CB4507')] <- 'N2'
          y$wild = y$haplotype!='N2'
          blens = do.call(rbind, mclapply(split(y, y$line), mc.cores=np, function(z) {rlei = rle(z$haplotype); rlew = rle(z$wild); data.frame(line = z$line, posr = z$posr, meanHapLen = rep(rlei$lengths/adjUnit, times = rlei$lengths), wild = rep(rlew$lengths/adjUnit, times = rlew$lengths))}))
          blens$wild = blens$meanHapLen!=blens$wild
          cbind(chrom=y$chrom[1], aggregate(data = blens, meanHapLen ~ posr+wild, mean))
        })))
}))
hgmap <- aggregate(data = haplos, pos ~ posr + chrom, function(x) mean(range(x)))
mhlw <- merge(mhlw, hgmap)

# plotting fn
plotROIhaplotypes <- function(roiChrom, roi, fnd, rilhaplotypes, hapLens, opref, buff=5e6, max_cM=10){
  
  # buff/max_cM = buffer in physical (bp)/genetic (cM) units around each gene
  
  ix = grep('AB1', names(fnd))
  fnd = cbind(fnd[,c('chrom', 'pos', 'genetic')], fnd[,ix:(ix+15)])
  pref = sprintf('%s_chr%s_%s', opref, roiChrom, roi)
  
  # founder haplotypes
  fndroi <- subset(fnd, chrom==roiChrom & pos >(roi-buff) & pos <(roi+buff))
  while(diff(range(fndroi$genetic)) > max_cM){
    buff = buff-1e4
    fndroi <- subset(fnd, chrom==roiChrom & pos >(roi-buff) & pos <(roi+buff))
  }
  
  hk <- scale(t(fndroi[,-(1:3)]))
  hk = tcrossprod(hk)
  hc = hclust(cluster::daisy(hk, 'manhattan'), method='ward.D2') 
  fndroi <- melt(fndroi, 1:3, value.name = 'genotype')
  fndroi$ix <- as.numeric(as.character(factor(fndroi$pos, labels = 1:length(unique(fndroi$pos)))))
  fndroi$variable <- factor(fndroi$variable, levels = hc$labels[hc$order])
  roix = fndroi$ix[which.min(abs(fndroi$pos-roi))]
  rocm = fndroi$genetic[which.min(abs(fndroi$pos-roi))]
  
  p <- ggplot(fndroi, aes(ix, variable, fill = factor(genotype))) + geom_tile() + scale_fill_discrete('') + labs(x='SNP', y='') + geom_vline(xintercept = roix, alpha=0.2, col = 'white')
  if(!roi %in% fndroi$pos) {fndroi <- subset(fndroi, ix!=roix); roix = fndroi$ix[which.min(abs(fndroi$pos-roi))]; p <- p + geom_vline(xintercept = roix, alpha=0.2, col = 'white')}
  ggsave(sprintf('%s_founderHaplotypes.png', pref), h = 4, w = 6)
  
  # RIL haplotypes - split by population and focal genotype
  shaplos = subset(rilhaplotypes, chrom==roiChrom & pos>(roi-buff) & pos<(roi+buff))
  # merge N2/CB4507
  shaplos$haplotype <- as.character(shaplos$haplotype)
  shaplos$haplotype[shaplos$haplotype %in% c("N2anc", "CB4507")] <- "N2/CB4507"
  
  orderWithinPops <- function(haps){
    # order lines by haplotype similarity, within each pop
    haps$pop <- gsub('L', '', substr(haps$line, 1, 6))
    do.call(rbind, lapply(split(haps, haps$pop), function(x){
      hmat <- dcast(x[,c('line', 'haplotype', 'pos')], pos~line, value.var='haplotype')[,-1]
      cdist <- function(X){
        # generic % identity fn
        n <- ncol(X); m <- nrow(X); K <- diag(n)
        for (i in 2:n) { for (j in 1:(i-1)) {K[i,j] <- K[j,i] <- sum(X[,i]==X[,j], na.rm=T)/m}}
        K
      }
      hk <- cdist(hmat)
      colnames(hk) <- colnames(hmat)
      diag(hk) <- 0
      hc = hclust(cluster::daisy(hk, 'manhattan'), method='ward.D2') 
      x$ix <- factorToInt(factor(x$line, labels = order(hc$order)))
      x
    }))
  }
  
  shaplos <- orderWithinPops(shaplos)

  # founder freq ~ pop
  ffreq = shaplos %>% group_by(chrom, posr, haplotype, pop) %>% count()
  popn <- aggregate(data = ffreq, n~pop+posr, sum)
  names(popn)[3] = 'N'
  ffreq <- merge(ffreq, popn)
  ffreq$sh = ffreq$haplotype=='N2/CB4507'
  ffreq$pop <- factor(ffreq$pop)
  ffreq$pop <- factor(ffreq$pop, levels = unique(ffreq$pop)[c(1,3,5,7,2,4,6,8,9,10)])
  p <- ggplot(ffreq, aes(posr, n/N, col = haplotype)) + geom_point(aes(shape=sh), stroke=0) + geom_line() + theme_classic() +
    geom_vline(xintercept = rocm) + facet_grid(pop~., scales='free', space='free') +
    scale_color_manual(values = pal3) + guides(colour = guide_legend("", override.aes = list(size=6))) + labs(x = 'cM') + scale_shape_discrete(guide=F)
  ggsave(sprintf('%s_sumHaplotypes.png', pref), h = 10, w = 4)
  
  # haplotype frequency over time
  ffreq$g <- 0; ffreq$g[grep('A6', ffreq$pop, invert = T)] <- as.numeric(substr(ffreq$pop[grep('A6', ffreq$pop, invert = T)], 4, 6))
  traj <- do.call(rbind, lapply(split(ffreq, list(ffreq$posr, ffreq$haplotype)), function(x) {
    if(len(unique(x$g))>2){
      # limit to focal CAs
      x = x[grep('GA', x$pop, invert = T),]
      fit = glm(cbind(n,N)~g, x, family = 'binomial')
      cf = coef(summary(fit))[2,c(1,3)]
      cbind(x[1,c('posr', 'chrom', 'haplotype')], coef = cf[1], z = cf[2])
    }
  }))
  traj$sh <- traj$haplotype=='N2/CB4507'
  p <- ggplot(traj, aes(posr, z, col = haplotype)) + geom_point(aes(shape=sh), stroke=0) + geom_line() + theme_classic() +
    geom_vline(xintercept = rocm) + scale_color_manual(values = pal3) + guides(colour = guide_legend("", override.aes = list(size=6, shape=15))) + labs(x = 'cM') + scale_shape_discrete(guide=F)
  ggsave(sprintf('%s_trajHaplotypes.png', pref), h = 4, w = 5)
  
  labs <- ffreq; labs$haplotype[labs$haplotype!='N2/CB4507'] <- 'wild'
  labs$popr = substr(labs$pop, 1, 3)
  labs$rep = substr(labs$pop, 3, 3)
  labs$popg = substr(labs$pop, 1, 2)
  labs <- do.call(rbind, lapply(split(labs, list(labs$posr, labs$popr, labs$g)), function(x) {if(nrow(x)>0) cbind(x[1,c('chrom', 'posr', 'g', 'popr', 'popg', 'rep')], wild = sum(x$n[x$haplotype=='wild'])/sum(x$n))}))
  p <- ggplot(labs, aes(posr, wild, col = popg, group = popr, linetype=rep)) + geom_point(stroke=0) + geom_line() + theme_classic() + 
    geom_vline(xintercept = rocm) + labs(x = 'cM', y = 'Frequency of non-N2 haplotypes') + facet_grid(.~g) + scale_color_discrete("") + theme(legend.position = 'top') + ylim(c(0,1))
  ggsave(sprintf('%s_trajWildHaplotypes.png', pref), h = 4, w = 5)
  # focal A6140 AND CA50/100 populations only
  p <- ggplot(subset(labs, popg!='GA'), aes(posr, wild, group = popr, linetype=rep)) + geom_point(stroke=0, alpha=0.3) + geom_line() + 
    theme_classic() + geom_vline(xintercept = rocm) + labs(x = 'cM', y = 'Frequency of non-N2 haplotypes') + facet_grid(.~g) + scale_color_discrete("") + theme(legend.position = 'top') + ylim(c(0,1)) + scale_x_continuous(breaks = scales::pretty_breaks(n=2))
  ggsave(sprintf('%s_trajWildHaplotypes-GA.png', pref), h = 4, w = 5)
  
  # mean haplotype length
  hroi <- subset(hapLens, chrom==roiChrom & pos >(roi-buff) & pos <(roi+buff))
  hroi$genetic = hroi$posr
  hroi$popr = substr(hroi$pop, 1, 3)
  hroi$rep = substr(hroi$pop, 3, 3)
  hroi$popg = substr(hroi$pop, 1, 2)
  hroi$g <- 0; hroi$g[grep('A6', hroi$pop, invert = T)] <- as.numeric(substr(hroi$pop[grep('A6', hroi$pop, invert = T)], 4, 6))
  hroi$wild = factor(hroi$wild, labels = c('N2', 'wild'))
  while(diff(range(hroi$genetic)) > max_cM){
    buff = buff-1e4
    hroi <- subset(hroi, pos >(roi-buff) & pos <(roi+buff))
  }
  p <- ggplot(hroi, aes(posr, meanHapLen, col = popg, group = popr, linetype=rep)) + geom_point(stroke=0) + geom_line() + 
    theme_classic() + geom_vline(xintercept = rocm) + labs(x = 'cM', y = 'Haplotype length (cM)') + facet_grid(wild~g) + scale_color_discrete("") + theme(legend.position = 'top')
  ggsave(sprintf('%s_hlenWildHaplotypes.png', pref), h = 4, w = 5)
  
}

sapply(1:nrow(goi), function(i) plotROIhaplotypes(roiChrom = goi$chrom[i], roi = goi$mid[i], fnd, rilhaplotypes = haplos, hapLens = mhlw, opref = goi$gene[i]))

