#!/usr/bin/env Rscript

#######################################################################
# LMN Aug 2019
# conversion of parsed Choreography output in R data container to WCON

# time, x/y, angular(turningRate)/area/bias/curve(curvature)/length/speed(velocity)/width (where present, in @MWT slot)
# see https://github.com/openworm/tracker-commons

# pass a config file (delimited, headered) listing per line:
# (1) a listing of .RData files (Parsed output from Choreography)
# (2) a prefix for the .wcon.zip (if blank/NA, the input file trimmed of prefix/suffix is used)
# (3:N) metadata (with names taken from the header).
# For CeMEE: 
#   sample      RIL/population ID
#   sample_type RIL, population, founder, MA etc
#   who         experimenter
#   place       Paris/Lisbon
#   date        YYYYMMDD
#   time        HHMMSS
#   media       NGM/NaCl
#######################################################################

require(jsonlite)
require(data.table)
require(parallel)
require(lubridate)

args = commandArgs(trailingOnly = T)
cfg  = args[1]  # config file
outd = args[2]  # all output written to common dir 

# unarray unit values recursively
unb <- function(i){
  if(length(i)==1){
    unbox(i)
  } else {
    lapply(i, function(j) unb(j))
  }
}

# common metadata tags
metal <- list(
  lab=list(
    name='EEV', 
    location='Institut de Biologie de l’École Normale Supérieure, CNRS UMR 8197, Inserm U1024, PSL 7 Research University, F-75005 Paris, France'
  ),
  arena=list(
    style='petri',
    size=90
  ),
  stage='adult'
)

# common units list for metadata and data
unitl = lapply(
  list(t='s',
       x='mm',
       y='mm',
       length='mm',
       width='mm',
       area='mm^2',
       speed='mm/s',
       curve='r',
       angular='r/s',
       bias='',
       food='HT115',
       software=list(name='MWT', 
                     version='1.3.0_r1035',
                     settings='-S --shadowless -q --plugin Reoutline --plugin Respine -N all',
                     featureID='@MWT'
       )
  ), 
  function(i) unb(i)
)

extractRda <- function(rda){
    
  load(rda)
  res = results[[3]]
  names(res)[names(res)=='XPosition'] <- 'x'
  names(res)[names(res)=='Yposition'] <- 'y'
  res <- na.exclude(res[,c('Time', 'IndWorm', 'x', 'y', 'Velocity', 'Bias', 'Area', 'Length', 'Width', 'Curvature', 'TurningRate')])
  
  lapply(split(res, res$IndWorm), function(b) {
    
    list(id = unbox(b$IndWorm[1]),
         t = list(b$Time),
         x = list(b$x),
         y = list(b$y),
         "@MWT" = list(angular = list(b$TurningRate),
                       area = list(b$Area),
                       bias = list(b$Bias),
                       curve = list(b$Curvature),
                       length = list(b$Length),
                       speed = list(b$Velocity),
                       width = list(b$Width)
         )
    )
    
  })
  
}

main <- function(cfg, outd, np=1){
  
  setwd(outd)
  config = fread(cfg, data.table = F, colClasses = c(date='character', time='character'))
  nda = nrow(config)
  
  config = config[order(config$sample_type, config$sample),]
  rdas = config[,1]
  prefs = config[,2]
  metadata = config[,3:ncol(config)]
  tstamp = ymd_hms(paste(metadata$date,metadata$time))
  tstamp = sapply(1:nda, function(i) strftime(tstamp[i], '%Y-%m-%dT%H:%M:%SZ', tz = sprintf('Europe/%s', metadata$place[i])))
  metadata$timestamp = tstamp
  metadata = metadata[,!names(metadata) %in% c('date', 'time')]
  
  o = mclapply(1:nda, mc.cores = np, function(i) {
    
    rda = rdas[i]
    pref = ifelse(is.na(prefs[i]), gsub('Parsed_data-', '', gsub('.RData', '', basename(rda))), prefs[i])
    outf = sprintf('%s.wcon', pref)
    osha = sprintf('%s.sha', outf)
    ozip = sprintf('%s.zip', outf)
    
    if(!file.exists(ozip)){
      cat(sprintf('converting %s\n', rda))
      
      oi = list(
        metadata = unb(c(metal, metadata[i,], file = rda)),
        units = unitl,
        data = extractRda(f)
      )
      
      write_json(oi, outf, pretty=T)
      system(sprintf('shasum %s > %s', outf, osha))
      zip(ozip, files = c(outf, osha))
      system(sprintf('rm %s %s', outf, osha))  
    } else {
      cat(sprintf('skipping %s\n', rda))
    }
    
  })
  
}

cfg = '~/rdas.cfg'
outd = '~/MWTdata/wcon'
main(cfg, outd, np=20)


makeRdaCfg <- function(){
  
  # thor
  # make cfg for all parsed data
  
  pdir = system(sprintf('find %s -name "Parsed*RData"', '/users/gev/mallard/Data/Cemee_data/'), intern=T)
  pdir = c(pdir, system('find /users/gev/mallard/Data/Parsed_behavior/B[0-9]* -name "Parsed*.RData"', intern = T))
  
  rdas = data.frame(rda=pdir, stringsAsFactors = F)
  bn = basename(rdas$rda)
  sp = tstrsplit(bn, '-')
  dd = gsub('_proc', '', sp[[3]])
  rdas$date = tstrsplit(dd, '_')[[1]]
  rdas$time = substr(tstrsplit(dd, '_')[[2]], 1, 6)
  rdas$block = sp[[2]]
  mm = gsub('.RData', '', sp[[4]])
  spp = tstrsplit(mm, '_')
  rdas$id = spp[[1]]
  rdas$media = spp[[2]]
  rdas = subset(rdas, media!='control')
  rdas = subset(rdas, block!='MaleControlBlock')
  rdas$media[rdas$media=='NGMl'] = 'NGM'
  # experimenter
  rdas$year = as.numeric(substr(rdas$date, 1, 4))
  rdas <- subset(rdas, year!=2015)
  rdas$who = 'BA'
  rdas$who[year==2016] <- 'AC'
  rdas$who[year %in% c(2018, 2019)] <- 'FM'
  # place
  rdas$place = 'Paris'
  rdas$place[rdas$year == 2012] = 'Lisbon'
  # sample/type
  rdas$sample = rdas$id = gsub('L0$', '', rdas$id)
  rdas$sample_type = NA
  rdas$sample[rdas$id=='N2anc' & !rdas$block %in% c('B302','B305')] = 'N2'
  # reformat
  x = "G50MOVA"; ix = grep(x, rdas$sample)
  rdas$sample[ix] <- gsub(x, '', rdas$sample[ix])
  rdas$sample[ix] <- sprintf('GA%s50%s', substr(rdas$sample[ix], 1, 1), substr(rdas$sample[ix], 2, 5))
  x = "G140A6"; ix = grep(x, rdas$sample)
  rdas$sample[ix] <- sprintf('A6140%s', gsub(x, '', rdas$sample[ix]))
  rdas$sample_type[grep('^PB250', rdas$sample)] <- 'MA line'
  rdas$sample_type[grep('^N2250', rdas$sample)] <- 'MA line'
  rdas$sample_type[rdas$sample=='N2'] = 'MA founder'
  ix = grep('L', rdas$sample, invert = T); ix = ix[is.na(rdas$sample_type[ix])]
  rdas$sample_type[ix] = 'population'
  cemeeFounders <- sort(c('AB1', 'CB4507', 'CB4855', 'CB4856', 'CB4858', 'MY1', 'MY16', 'JU400', 'JU319', 'RC301', 'PX179', 'N2anc', 'PB306', 'PX174', 'CB4852', 'JU345'))
  rdas$sample_type[rdas$sample %in% cemeeFounders] <- 'founder'
  rdas$sample[rdas$id=='PBanc'] = 'PB306'
  rdas$sample_type[rdas$id=='PBanc' & rdas$sample=='PB306'] = 'MA founder'
  rdas$sample_type[is.na(rdas$sample_type)] = 'RIL'
  # Patrick's lines
  rdas$sample_type[rdas$sample == 'CX12311'] = 'NIL'
  rdas$sample_type[rdas$sample == 'LSJ2'] = 'mutant'
  # exc. male depleted samples
  rdas = rdas[grep('noM', rdas$sample, invert = T),]
  # Bruno's old N2
  rdas <- subset(rdas, !(who=='BA' & sample=='N2'))
  
  rdas$pref=sprintf("%s_%s_%s", rdas$sample, rdas$date, rdas$time)
  rdas = rdas[,c('rda','pref','sample','sample_type','who','place','date','time','media')]
  rdas = rdas[order(as.numeric(rdas$date), as.numeric(rdas$time)),]
  
  # there is data split by sex for just one block (keep original)
  rdas = rdas[grep('splitted', rdas$rda, invert=T),]
  
  fwrite(rdas, '~/rdas.cfg')
  
  # names(bt)[names(bt) %in% metad$data_group_name]
  # names(bt)[!names(bt) %in% metad$data_group_name]
  # meh ignore metadata. get line/food for blobs filenames
  # metad <- fread('~/Documents/cemee/MWT2WCON/Luke_header.txt')
  # names(metad)[6:7] = c('date', 'time')
  
}

