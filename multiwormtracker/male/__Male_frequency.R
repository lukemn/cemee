#!/usr/bin/env Rscript
#
# June 2019, LMN
#
# Classify parsed Choreography tracks using an extreme gradient boosting model
# trained on populations (including mono & dioecious) and inbred founders (extreme 40 of 70 samples, -CB4856).
# The current model is saved in xgb_preds.rda.
# Track-level locomotion statistics are calculated (locostat), 
# sex is predicted and annotated tracks are dumped to .rda (assignSexToTracks), 
# then male/herm frequencies are sampled from small time slices (sampleHMs).

require(xgboost)
require(parallel)
require(data.table)

args  = commandArgs(trailingOnly=T)
WD    = args[1] # Working directory containing parsed Choreography files (globbed as `Parsed*.RData`). 
                # Files for tracks passing QC will be saved in `WD/sex/`, with a new column `sex`.
PREF  = args[2] # Prefix for saving intermediate summary traits (`WD/sex/PREF_simpleTraits.RData`)
XGBmod= args[3] # trained XBG model for sex prediction
FNP   = 8       # Default parallel threads across/within files
NP    = 2

## Sanity checks.
# WD must exist and be writable
stopifnot(file.exists(WD))
stopifnot(system(sprintf('mkdir %s/sex', WD), ignore.stderr = F)==0)
# adjust threads if not running on Thor
NC = detectCores()
if(NC==1) FNP=NP=1
while(FNP*NP > NC) FNP=FNP-1
if(file.exists(XGBmod)){
  load(XGBmod)
} else {
  cat(sprintf("Can't find classification model %s\n", XGBmod))
  stop()
}

locostat <- function(D, minExpDuration=10, Lcensor=5, subsampleFrameRate=4, minTracks=30, minTrackDurationSec=3, minTrackObsPerState=3, agSamples = c(10, 15, 20), NP=1){
  
  # Generate simple summary traits from MWT Choreography output.
  # Returns NULL if sample fails QC.
  #
  # Parameters, in order of application:
  # minExpDuration (minutes)
  # Lcensor: discard early data (minutes)
  # subsampleFrameRate: subsample raw track data (Hertz)
  # minTracks: minimum number of unique tracks to analyse
  # minTrackDurationSec: subsampled tracks are filtered on total length (seconds)
  # minTrackObsPerState: runs of consecutive Bias states within tracks are filtered on length (Fwd|Back >1), and length of longest run (Fwd|Back >=minTrackObsPerState)
  # agSamples: times (minutes) at which data are tested for non-random worm aggregation within the imaging field (see agQ() below)
  
  load(D, verbose=T)
  res <- results[[3]]
  names(res)[names(res)=='XPosition'] <- 'Xposition'
  res <- res[,c('Time', 'IndWorm', 'Xposition', 'Yposition', 'Velocity', 'Bias', 'Persistence', 'Area', 'Length', 'Width', 'Curvature')]
  
  if(max(res$Time) >= minExpDuration*60){
    
    N_raw = length(unique(res$IndWorm))
    # discard first Lcensor minutes
    res <- subset(res, Time > Lcensor*60)
    # order (should be unneccesary)
    res <- res[order(res$IndWorm, res$Time),]
    # subsample to mean subsampleFrameRate Hz / track
    res$w <- round(res$Time / (10/subsampleFrameRate),1) * (10/subsampleFrameRate)
    res <- do.call(rbind, mclapply(split(res, res$IndWorm), mc.cores = NP, function(x) x[!duplicated(cbind(x$w, x$Bias)),]))
    rownames(res) <- NULL
    res <- na.exclude(res)
    tracks = unique(res[,c('IndWorm', 'Persistence')])
    Ntracks = nrow(tracks)
    # approximations to plate density
    logd = log(sum(tracks$Persistence)/max(res$Time))
    # mean tracks / s
    trackPerS <- do.call(rbind, mclapply(split(res, res$IndWorm), mc.cores = NP, function(x) x[!duplicated(round(x$Time)),]))
    logmtps = log(mean(as.numeric(table(round(trackPerS$Time)))))
    
    agQ <- function(df, t, xl, yl, NP, nsamp=5e3, fgrid=1e3){
      # empirical quantiles for mean pairwise distance between all worms at time t minutes
      # based on nsamp samples across a field grid defined by the positions of extreme worms
      lt <- function(x) {x = as.matrix(x); x[lower.tri(x, diag = F)]}
      dsamp <- subset(df, w==t*60)
      if(nrow(dsamp) >= minTracks){
        dsamp <- dsamp[!duplicated(dsamp$IndWorm),c('Xposition', 'Yposition')]
        dm <- as.matrix(dist(dsamp))
        ag_obs = mean(lt(dm))
        n = nrow(dm)
        ag_exp = unlist(mclapply(1:nsamp, mc.cores = NP, function(i) {
          si = cbind(sample(seq(xl[1], xl[2], length.out = fgrid), n, replace=T), sample(seq(yl[1], yl[2], length.out = fgrid), n, replace=T))
          mean(lt(dist(si)))
        }))
        sum(ag_obs > ag_exp)/nsamp 
      }
    }

    if(Ntracks >= minTracks){
      # observed imaging field
      xl = range(res$Xposition)
      yl = range(res$Yposition)  
      mmedian <- function(x) ifelse(length(x)>1, median(x, na.rm=T), NA)
      mvar <- function(x) ifelse(length(x)>1, var(x, na.rm=T), NA)
      tstat = do.call(rbind, mclapply(split(res, res$IndWorm), mc.cores=NP, function(x) {
        f = subset(x, Bias==1)
        s = subset(x, Bias==0)
        b = subset(x, Bias== -1)
        r = rle(x$Bias)
        xr = split(x, rep.int(1:length(r$values), times = r$lengths))
        # filter subsampled tracks on total time (>minTrackDurationSec)
        # filter runs within tracks on length (f|b >1), length of longest run (f|b >=minTrackObsPerState)
        # variances may be NA
        xrf = xr[r$values==1]
        xrb = xr[r$values== -1]
        xrfl = ifelse(length(xrf)>0, max(sapply(xrf, nrow)), 0)
        xrbl = ifelse(length(xrb)>0, max(sapply(xrb, nrow)), 0)
        
        if(all(diff(range(x$Time))>minTrackDurationSec & (length(xrf)>1 | length(xrb)>1) & (xrfl>=minTrackObsPerState | xrbl>=minTrackObsPerState))){
          data.frame(area.F = mmedian(f$Area),
                     area.F.var = mvar(f$Area),
                     area.S = mmedian(s$Area),
                     area.S.var = mvar(s$Area),
                     length.F = mmedian(f$Length),
                     width.F = mmedian(f$Width),
                     length.F.var = mvar(f$Length),
                     width.F.var = mvar(f$Width),
                     length.S = mmedian(s$Length),
                     width.S = mmedian(s$Width),
                     length.S.var = mvar(s$Length),
                     width.S.var = mvar(s$Width),
                     velocity = mmedian(x$Velocity),
                     velocity.var = mvar(x$Velocity),
                     velocity.F = mmedian(f$Velocity),
                     velocity.B = mmedian(b$Velocity),
                     velocity.F.var = mvar(f$Velocity),
                     velocity.B.var = mvar(b$Velocity),               
                     acceleration = mmedian(diff(x$Velocity)),
                     acceleration.var = mvar(diff(x$Velocity)),
                     acceleration.F = mmedian(unlist(lapply(xrf, function(y) diff(y$Velocity)))),
                     acceleration.B = mmedian(unlist(lapply(xrb, function(y) diff(y$Velocity)))),
                     acceleration.F.var = mvar(unlist(lapply(xrf, function(y) diff(y$Velocity)))),
                     acceleration.B.var = mvar(unlist(lapply(xrb, function(y) diff(y$Velocity)))),
                     run = mmedian(r$lengths),
                     run.var = mvar(r$lengths),
                     run.F = mmedian(unlist(lapply(xrf, nrow))),
                     run.B = mmedian(unlist(lapply(xrb, nrow))),
                     run.F.var = mvar(unlist(lapply(xrf, nrow))),
                     run.B.var = mvar(unlist(lapply(xrb, nrow))),
                     curvature = mmedian(x$Curvature),
                     curvature.var = mvar(x$Curvature),
                     curvature.F = mmedian(f$Curvature),
                     curvature.B = mmedian(b$Curvature),
                     curvature.S = mmedian(s$Curvature),
                     curvature.F.var = mvar(f$Curvature),
                     curvature.B.var = mvar(b$Curvature),
                     curvature.S.var = mvar(s$Curvature),
                     persistence = x$Persistence[1],
                     worm = x$IndWorm[1],
                     trackStartTime = x$Time[1]
          ) 
        }
      }))
      
      agqs <- sapply(agSamples, function(t) agQ(res, t, xl, yl, NP))
      names(agqs) = paste0('aggregation.t', agSamples)
      
      list(trackStats = tstat, 
           aggregationQ = agqs, 
           N_tracks_raw = N_raw,
           N_tracks_pass = nrow(tstat),
           logDensity = logd,
           logMeanTracksPerSecond = logmtps,
           L_pass = diff(range(res$Time)),
           meanObsPerTrack_pass = mean(as.numeric(table(res$IndWorm))),
           R_velocityTime = cor(res$Velocity, res$Time),
           file = D
      ) 
    }
  }
}

assignSexToTracks <- function(dfo, parsed_files){

  # From processed track-level, summary traits for all samples (data.frame dfo) generated by locostat,
  # classify the Choreography tracks using the xbg model.
  
  sexid = xpreds
  sexid$sex = xpreds$psex

  # enforce consistent log transformations on raw traits between classifier and new data
  # this is the matrix of track-level simple traits for all samples
  xtraits = xmod$feature_names
  Xd = dfo
  for(i in names(Xd)) if(paste0('ln.', i) %in% xtraits) {Xd[,i] = log(Xd[,i]); names(Xd)[names(Xd)==i] <- paste0('ln.', i)}
  stopifnot(sum(names(Xd) %in% xtraits)==length(xtraits))
  Xd = as.matrix(Xd[,xtraits])
  dfo$sex = predict(xmod, Xd)
  # binarize probabilities
  dfo$sex = factor(dfo$sex>0.5, labels = c('male', 'herm'))

  # return to parsed tracks (one per sample), merge and dump.
  # merging is based on formatted filename.
  for(i in parsed_files){
    xidi = sapply(strsplit(basename(i), '-'), function(x) paste(x[2:3], collapse='-'))
    if(xidi %in% dfo$xid){
      load(i, verbose=T)
      res <- results[[3]]
      subx <- subset(sexid, sexid$xid==xidi)[,c('worm', 'sex')]
      names(subx) <- c('IndWorm', 'sex')
      reso <- merge(res, subx, sort=F)
      results[[3]] <- reso[order(reso$IndWorm, reso$Time),]
      of = sprintf('%s/sex/%s', WD, gsub('.RData', '_classified.RData', basename(i)))
      save(results, file = of)
    }
  } 
}

sampleHMs <- function(classifiedTrackFiles, sampleStartMin=5, sampleEndMin=20, sampleFreqPerMin=2, sampleWindowS=1, offset=0, subsampleFrameRate=4, minTrackLengthS=10, NP=12){
  
  # Sample herm/male frequencies from sex-assigned tracks to account for track frequency/sex confounding.
  #
  # Tracks are subsampled to 4Hz, filtered to minTrackLengthS seconds (~.1% quantile), 
  # and dumped to a single (potentially very large) file.
  # Samples of sampleWindowS (seconds) are taken at sampleFreqPerMin over sampleStartMin:sampleEndMin (+offset) intervals,
  # (potential conflict with locostat Lcensor parameter if modified from the defaults).
  # output is a data.frame of frequencies per time slice and experimental file (defined by date/time of acquisition), 
  # ready for estimation by glm across replicates.
  
  # sampleStartMin=2; sampleEndMin=20; sampleFreqPerMin=2; sampleWindowS=1; offset=0; subsampleFrameRate=4; minTrackLengthS=10; NP=12
  
  # filter and subsample tracks
  o = sprintf('%s/sex/%s_mergedClassifiedTracks.rda', WD, PREF)
  if(!file.exists(o)){
    tracks <- do.call(rbind, lapply(classifiedTrackFiles, function(f) {
      load(f, verbose=T)
      res = results[[3]]
      res <- res[order(res$IndWorm, res$Time),]
      res = subset(res, Persistence>minTrackLengthS)
      # viva la France
      res$Time = as.numeric(tstrsplit(res$Time, ',')[[1]])
      if(nrow(res)>0){
        res$w <- round(res$Time / (10/subsampleFrameRate),1) * (10/subsampleFrameRate)
        res <- do.call(rbind, mclapply(split(res, res$IndWorm), mc.cores = NP, function(x) x[!duplicated(x$w),]))
        res$id = results[[2]]
        handle = tstrsplit(basename(f), '_')
        res$date = tstrsplit(handle[[2]], '-')[[3]]
        res$xtime = handle[[3]]
        res  
      }
    }))
    tracks$sex <- factor(tracks$sex)
    save(tracks, file=o)
  } else {
    load(o)
  }
  
  # sample male/herm freqs for each MWT file. NB this eats RAM
  w = seq(sampleStartMin, sampleEndMin, 1/sampleFreqPerMin)*60
  sampledFreqs <- do.call(rbind, mclapply(split(tracks, paste(tracks$date, tracks$xtime)), mc.cores = 1, function(x) {
    mfreqs = do.call(rbind, lapply(w, function(i) table(subset(x, round((Time+offset)/sampleWindowS)==round(i/sampleWindowS))$sex)))
    cbind(x[1,c('id', 'xtime', 'date')], h=mfreqs[,1], m=mfreqs[,2], w=w[1:nrow(mfreqs)], row.names=NULL)
  }))
  save(sampledFreqs, file = sprintf('%s/sex/%s_sampledFreqs_w%s.rda', WD, PREF, sampleWindowS))
}

main <- function(){
  
  # extract simple traits from each set of parsed Choreography tracks
  pfiles = Sys.glob(sprintf("%s/Parsed*.RData", WD))
  out = mclapply(pfiles, mc.cores = FNP, function(f) tryCatch(locostat(f, NP=NP), error = function(e) f))
  outf = sprintf('%s/%s_simpleTraits.RData', WD, PREF)
  save(out, pfiles, file = outf)
  
  # merge tracks that pass QC
  nulls <- unlist(lapply(seq_along(out), function(x) is.null(out[[x]][[1]])))
  errs <- unlist(lapply(seq_along(out), function(x) length(out[[x]])))==1
  print(sprintf('%s: %s (of %s) failed QC, %s errors', basename(outf), sum(nulls), length(nulls), sum(errs)))
  out <- out[!(nulls|errs)]
  pfiles <- pfiles[!(nulls|errs)]
  ids = unlist(lapply(strsplit(basename(pfiles), '-'), function(x) paste(x[2:3], collapse='-')))
  merged <- do.call(rbind, lapply(seq_along(ids), function(x) cbind(xid = ids[x], out[[x]][[1]], out[[x]][3:6])))
  
  # predict sex for each track and dump per sample .RData in WD/sex
  assignSexToTracks(merged, pfiles)
  
  # sample male/herm frequencies from the classified tracks and dump single .rda in WD/sex
  cfiles = Sys.glob(sprintf("%s/sex/Parsed*.RData", WD))
  sampleHMs(cfiles, NP=FNP*NP)
}





