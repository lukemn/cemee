#!/usr/bin/env Rscript

#######################################################################
# LMN Aug 2019
# conversion of MWT blobs to WCON
# https://github.com/openworm/tracker-commons/blob/master/WCON_format.md

# time, x/y, length/width (in @MWT slot) at minimum
# by default, filtered to only tracks with data for (11-point) skeleton and outline
# pixels converted to mm, assuming 0.026mm per pixel
# outlines (walks) are unaltered from Rex's compressed pixel-based form:
# walk = [{px=[x & y of starting pixel, pixel length=1], n=number of pixels in outline, '4'=base64 encoded binary array giving x/y sign from starting pixel}]

# pass a config file (delimited, headered) listing per line:
# (1) a directory containing blobs files,
# (2) a prefix for the .wcon.zip output (if blank/NA, the parent directory is used)
# (3:N) metadata (with names taken from the header).
# For CeMEE: 
#   sample      unique ID
#   sample_type RIL, NIL, population, founder, MA line, MA founder etc
#   who         experimenter
#   place       Paris/Lisbon
#   date        YYYYMMDD
#   time        HHMMSS
#   media       NGM/NaCl

# output is written to the each blobs directory, then zipped along with a sha1 hash
###################################################################################

require(jsonlite)
require(data.table)
require(parallel)
require(lubridate)
require(dplyr)

args = commandArgs(trailingOnly = T)
cfg  = args[1] # config file 

# pixel > mm conversion
mmperpixel = 0.026
# unarray unit values recursively
unb <- function(i){
  if(length(i)==1){
    unbox(i)
  } else {
    lapply(i, function(j) unb(j))
  }
}

# common metadata tags
metal <- list(lab="EEV")

# common units list for metadata and data
unitl = lapply(
  list(t='s',
       x='mm',
       y='mm',
       length='mm',
       width='mm',
       food='HT115',
       software=list(name='MWT', version='1.3.0_r1035')), 
  function(i) unb(i)
)

parseblobs <- function(blobf, filterOnOutline=T){
  
  blockparser <- function(b, wormid){
    
    cols = lapply(b, function(l) {l=strsplit(l, ' ')[[1]]; l[l!='']})
    xes = lapply(cols, function(l) as.numeric(l[3])*mmperpixel)
    yes = lapply(cols, function(l) as.numeric(l[4])*mmperpixel)
    
    if(filterOnOutline){
    
      skeloffsets = lapply(cols, function(l) as.numeric(l[12:33])*mmperpixel)
      nl = length(b)
      
      list(id = unbox(wormid), 
           t = sapply(cols, function(l) as.numeric(l[2])),
           x = lapply(1:nl, function(i) xes[[i]] + skeloffsets[[i]][seq(1, 22, 2)]),
           y = lapply(1:nl, function(i) yes[[i]] + skeloffsets[[i]][seq(2, 22, 2)]),
           "@MWT" = list(length = sapply(cols, function(l) as.numeric(l[9])*mmperpixel),
                         width = sapply(cols, function(l) as.numeric(l[10])*mmperpixel)
           ),
           walk = list(px = lapply(cols, function(l) c(as.numeric(l[35:36]),1)),
                       n = sapply(cols, function(l) as.numeric(l[37])),
                       "4" = sapply(cols, function(l) l[38])
           )
      )
    } else {
      # just x,y,length,width
      list(id = unbox(wormid), 
           t = sapply(cols, function(l) as.numeric(l[2])),
           x = unlist(xes),
           y = unlist(yes),
           "@MWT" = list(length = sapply(cols, function(l) as.numeric(l[9])*mmperpixel),
                         width = sapply(cols, function(l) as.numeric(l[10])*mmperpixel)
           )
      )
    }
  }
  
  lines_ = readLines(blobf)
  ix = grep('^%', lines_)
  wormids = tstrsplit(lines_[ix], ' ')[[2]]
  
  parsed = lapply(seq_along(ix), function(i) {
    if(i==length(ix)){
      b = lines_[(ix[i]+1):length(lines_)]
    } else {
      b = lines_[(ix[i]+1):(ix[i+1]-1)]
    }
    tryCatch(blockparser(b, wormids[i]), error = function(e){NULL})  
  })
  
  if(filterOnOutline){
    # many lines do not have skeleton/outlines, ignoring them
    # lls = unlist(lapply(lines_, function(i) length(strsplit(i, ' ')[[1]])))
    # table(lls)
    isnull = sapply(parsed, function(i) is.null(unlist(i)))
    cat(sprintf('\t%s/%s tracks without skeleton/outline removed\n', sum(isnull), length(isnull)))
    data_ = parsed[!isnull]  
  } else {
    data_ = parsed
  }
  
  return(data_)
  
}

main <- function(cfg, np=1){
  
  config = fread(cfg, data.table = F, colClasses = c(date='character', time='character'))
  ndir = nrow(config)
  blobdir = config[,1]
  prefs = config[,2]
  metadata = config[,3:ncol(config)]
  tstamp = ymd_hms(paste(metadata$date,metadata$time))
  tstamp = sapply(1:ndir, function(i) strftime(tstamp[i], '%Y-%m-%dT%H:%M:%SZ', tz = sprintf('Europe/%s', metadata$place[i])))
  metadata$timestamp = tstamp
  metadata = metadata[,!names(metadata) %in% c('date', 'time')]
  
  o = lapply(1:ndir, function(i) {
    
    diri = blobdir[i]
    blobsf = Sys.glob(sprintf('%s/*.blobs', diri))
    pref = ifelse(is.na(prefs[i]), basename(diri), prefs[i])
    outf = sprintf('%s/%s.wcon', blobdir[i], pref)
    osha = sprintf('%s.sha', outf)
    ozip = sprintf('%s.zip', outf)
    
    if(!file.exists(ozip)){
      cat(sprintf('converting %s: %s blobs\n', diri, length(blobsf)))
      
      oi = list(
        metadata = unb(c(metal, metadata[i,])),
        units = unitl,
        data = unlist(mclapply(blobsf, mc.cores=np, function(f) parseblobs(f, filterOnOutline = F)), recursive = F)
      )
      
      write_json(oi, outf, pretty=T)
      system(sprintf('shasum %s > %s', outf, osha))
      zip(ozip, files = c(outf, osha))
      system(sprintf('rm %s %s', outf, osha))  
    } else {
      cat(sprintf('skipping %s\n', diri))
    }
    
  })
  
}

# cfg = '~/Documents/cemee/MWT2WCON/mwt_config'
cfg = '~/Documents/cemee/MWT2WCON/fm_blob.cfg'
main(cfg, np=1)

makeFMcfg <- function(){
  
  # make cfg (for the raw data that we have)
  blobsParentDir = '/Volumes/GEV_SAV_01/mallard/behavioral_data/Raw_behavior/'
  dirs = data.frame(dir = system(sprintf('find %s -type d -mindepth 2', blobsParentDir), intern=T), stringsAsFactors = F)
  # dirs <- fread('~/Documents/cemee/MWT2WCON/fm_blob_dirs', header = F); names(dirs) = 'dir'
  dirs = dirs[grep('^2', basename(dirs$dir)),,drop=F]
  dirs$date = tstrsplit(basename(dirs$dir), '_')[[1]]
  dirs$time = tstrsplit(basename(dirs$dir), '_')[[2]]
  dirs$block = basename(dirname(dirs$dir))
  (bt = table(dirs$block))
  # names(bt)[names(bt) %in% metad$data_group_name]
  # names(bt)[!names(bt) %in% metad$data_group_name]
  
  # meh ignore metadata. get line/food for blobs filenames
  # metad <- fread('~/Documents/cemee/MWT2WCON/Luke_header.txt')
  # names(metad)[6:7] = c('date', 'time')
  
  dirm = bind_rows(lapply(dirs$dir, function(i) {
    print(i)
    fs = basename(system(sprintf('ls %s/*blobs', i), intern = T, ignore.stderr = T))
    if(len(fs)>0){
      sp = tstrsplit(fs, '_')
      id = unique(sp[[1]])
      media = unique(sp[[2]])
      block = substr(unique(sp[[3]]), 1, 4)
      stopifnot(all(len(line)==1, len(food)==1, len(block)==1))
      data.frame(dir = i, id, media, block, stringsAsFactors = F)  
    }
  }))
  
  mdir <- merge(dirs, dirm)
  mdir$pref = ''
  mdir$who = 'FM'
  mdir$place = 'Paris'
  mdir$sample = gsub('L0', '', mdir$id)
  mdir$sample[mdir$id=='N2anc' & !mdir$block %in% c('B302','B305')] = 'N2'
  mdir$sample_type = NA
  mdir$sample_type[grep('^PB250', mdir$sample)] <- 'MA line'
  mdir$sample_type[grep('^N2250', mdir$sample)] <- 'MA line'
  mdir$sample_type[mdir$sample=='PBanc'] = 'MA founder'
  mdir$sample_type[mdir$sample=='N2'] = 'MA founder'
  ix = grep('L', mdir$sample, invert = T); ix = ix[is.na(mdir$sample_type[ix])]
  mdir$sample_type[ix] = 'population'
  mdir$sample_type[mdir$sample %in% cemeeFounders] <- 'founder'
  mdir$sample_type[is.na(mdir$sample_type)] = 'RIL'
  # exc. male depleted samples
  mdir = mdir[grep('noM', mdir$sample, invert = T),]
  mdir = mdir[,c('dir','pref','sample','sample_type','who', 'place','date','time','media')]
  fwrite(mdir, '~/Documents/cemee/MWT2WCON/fm_blob.cfg')
}

