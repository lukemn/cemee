#!/usr/bin/env Rscript

#######################################################################
# LMN Aug 2019
# conversion of MWT blobs to WCON

# length, width, 11-point skeleton, and outline retained
# pixels converted to mm, assuming 0.026mm per pixel
# outlines are unaltered from Rex's compressed pixel-based form and so are given as a length 4 string array:
# x & y of starting pixel, number of pixels in outline, then ASCII string giving bit pairs giving offset coord/direction + magnitudes from starting pixel)
# see https://github.com/openworm/tracker-commons
# {
#   "metadata":{
#     "who":"Open Worm",
#     "timestamp":"2016-01-22T17:44:48",
#     "protocol":"Numbers made up by hand for an example!"
#   },
#   "units":{"t":"s", "x":"mm", "y":"mm", "length":"mm", "width":"mm", "outline":"rex"},
#   "data":[
#     {"id":"worm1", "t":[0.1], "x":[[0.33, 0.65, 0.8, 1.1, 1.2, ...]], "y":[[2.31, 2.25, 2.0, 1.87, 1.66, ...]], "length":[[40.1]], "width":[[3.1]], "outline":[["11", "911", "118", "MMMeeMUGEEUUFFIIUVIVPc<<c<3<30`02<0228P0"]]},
#     {"id":"worm1", "t":[0.3], "x":[[0.27, 0.6, 0.75, 1.0, 1.1, ...]], "y":[[2.4, 2.3, 2.07, 1.78, 1.75, ...]]}, "length":[[64.3]], "width":[[13.2]], "outline":[["1040", "2303", "237", "MMGGMGEGEEeEMMEWedPX82080020202280P"]]}
#     ]
# }

# pass a config file (delimited, headered) listing per line:
# (1) a directory containing blobs files,
# (2) a prefix for the .wcon.zip (if blank/NA, the parent directory is used)
# (3:N) metadata (with names taken from the header).
# For CeMEE: 
#   sample      RIL/population ID
#   sample_type "RIL", "population"
#   date        YYYYMMDD
#   time        HHMMSS
#   who         experimenter
#   media       NGM/NaCl
#   temperature C
#   humidity    %
#######################################################################

require(jsonlite)
require(data.table)
require(parallel)
library(lubridate)

args = commandArgs(trailingOnly = T)
cfg  = args[1] # config file 

# pixel > mm conversion
mmperpixel = 0.026
# json converter
toj = jsonlite::toJSON
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
       outline='rex', 
       temperature='C',
       humidity='%',
       food='HT115',
       media='NGM',
       software=list(name='MWT', version='1.3.0_r1035')), 
  function(i) unb(i)
)

parseblobs <- function(blobf){
  
  blockparser <- function(b, wormid){
    
    nl = length(b)
    cols = lapply(b, function(l) {l=strsplit(l, ' ')[[1]]; l[l!='']})
    xes = lapply(cols, function(l) as.numeric(l[3])*mmperpixel)
    yes = lapply(cols, function(l) as.numeric(l[4])*mmperpixel)
    skeloffsets = lapply(cols, function(l) as.numeric(l[12:33])*mmperpixel)
    
    list(id = unbox(wormid), 
         t = sapply(cols, function(l) as.numeric(l[2])),
         x = lapply(1:nl, function(i) xes[[i]] + skeloffsets[[i]][seq(1, 22, 2)]),
         y = lapply(1:nl, function(i) yes[[i]] + skeloffsets[[i]][seq(2, 22, 2)]),
         length = sapply(cols, function(l) as.numeric(l[9])*mmperpixel),
         width = sapply(cols, function(l) as.numeric(l[10])*mmperpixel),
         outline = lapply(cols, function(l) c(as.numeric(l[35:37]), l[[38]]))
    )
  }
  
  lines_ = readLines(blobf)
  ix = grep('^%', lines_)
  wormids = tstrsplit(lines_[ix], ' ')[[2]]
  
  parsed = lapply(seq_along(ix), function(i) {
    if(i==length(ix)){
      datal = (ix[i]+1):length(lines_)
    } else {
      datal = (ix[i]+1):(ix[i+1]-1)
    }
    tryCatch(blockparser(lines_[datal], wormids[i]), error = function(e){NULL})  
  })
  
  # many lines do not have skeleton/outlines, ignoring them
  # lls = unlist(lapply(lines_, function(i) length(strsplit(i, ' ')[[1]])))
  # table(lls)
  isnull = sapply(parsed, function(i) is.null(unlist(i)))
  sprintf('%s incomplete data lines rmeoved')
  data_ = parsed[!isnull]
  return(data_)
}

main <- function(cfg, np=1){
  
  cfg = '~/Documents/cemee/MWT2WCON/mwt_config'
  config = fread(cfg, data.table = F)
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
    blobsf = Sys.glob(sprintf('%s/*.blobs', diri))[1:4]
    pref = ifelse(is.na(prefs[i]), basename(diri), prefs[i])
    outf = sprintf('%s/%s.wcon', blobdir[i], pref)
    
    oi = list(
      metadata = unb(c(metal, metadata[i,])),
      units = unitl,
      data = unlist(mclapply(blobsf, mc.cores=np, function(f) parseblobs(f)), recursive = F)
    )
    
    write_json(oi, outf, pretty=T)
    zip(paste0(outf, '.zip'), outf)
    system(sprintf('rm %s', outf))
  })
  
}

main(cfg, np=4)

