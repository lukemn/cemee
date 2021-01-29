


# list of simple traits
load('~/Downloads/Calibration3_simpleTraits.RData', verbose=T)
merged <- do.call(rbind, lapply(out, function(x) cbind(xid = basename(x$file[1]), x[[1]], x[3:6])))

xids = gsub('.RData', '', tstrsplit(merged$xid, '-')[[4]])
merged$env = substr(xids, 1, 1)
merged$ms = substr(xids, 3, 4)
merged$n = as.numeric(tstrsplit(xids, '_')[[3]])
merged$n[merged$n==100] = 1000

tracks = merged

XGB <- function(){}
require(xgboost)
require(mclust)

# first classify males/herms in dioceous using 3-traits 
subt=c('ln.area.F', 'ln.width.F', 'velocity.F')

# expect ~normal trait distributions
for(i in c('area.F', 'width.F')){
  tracks[,i] = log(tracks[,i])
  names(tracks)[names(tracks)==i] = paste0('ln.', names(tracks)[names(tracks)==i])
}

# remove NAs
tracks = subset(tracks, !(is.na(ln.area.F) | is.na(ln.width.F) | is.na(velocity.F)))

# split by mating system
nctrls = subset(tracks, ms %in% paste0('M', 1:3))
pctrls = subset(tracks, ms %in% paste0('D', 1:3))

ppreds <- do.call(rbind, mclapply(split(pctrls, pctrls$xid), mc.cores = 6, function(x) {
  m <- Mclust(x[,subt], G=1:2, modelNames = 'VII')
  # polarise by area. males are g1
  g = m$class; mg <- aggregate(data = data.frame(g, y=x[,'ln.area.F']), y~g, median); if(mg$y[1]>mg$y[2]) g = (g %% 2) + 1
  data.frame(x[,c(subt, c('xid', 'env', 'ms', 'n', 'worm'))], g, u=m$uncertainty)
}))

ggplot(ppreds, aes(u)) + geom_histogram(binwidth = .01) + geom_vline(aes(xintercept=quantile(u, 0.9)))

ppreds <- subset(ppreds, u < quantile(u, 0.90))
ppreds$sex <- factor(ppreds$g, labels = c('male', 'herm'))
dtrain <- rbind(ppreds[,c('xid', 'env', 'ms', 'n', 'sex', 'worm')], cbind(nctrls[,c('xid', 'env', 'ms', 'n', 'worm')], sex='herm'))
dtrain = merge(tracks, dtrain)
dtrain$env_n = as.numeric(as.character(factor(dtrain$env, labels=c(1, 1000))))
traits = names(tracks)[2:39]


# fit separate models for NGM, NaCl

dtrain = dtrain[dtrain$env=='N',]
dtrain = dtrain[dtrain$env=='S',]

Xd = xgb.DMatrix(as.matrix(dtrain[,traits]), label = factorToInt(factor(dtrain$sex, labels=c(0,1))))
param <- list(max_depth = 4, eta = 0.3, subsample = 0.8, silent = 1, nthread = 4, objective = "binary:logistic", eval_metric = "error")
# x validation
cv <- xgb.cv(data = Xd, params = param, nfold=3, nrounds=120)

xmod = xgboost(params = param, data = Xd, nrounds = 100, early_stopping_rounds = 10)
(fimp = xgb.importance(feature_names = traits, model=xmod))

xgb.save(xmod, fname = '/Users/noble/Documents/github/cemee/multiwormtracker/male/tp_ngm_xgb_preds')
xgb.save(xmod, fname = '/Users/noble/Documents/github/cemee/multiwormtracker/male/tp_nacl_xgb_preds')

# in sample predictions
th = 0.5
xpredt = predict(xmod, Xd)
table(dtrain$sex, xpredt>th)
sum(diag(table(dtrain$sex, xpredt>th)))/nrow(dtrain)

# out of sample
load('~/Downloads/EEhom_G5_simpleTraits.RData', verbose=T)
merged <- do.call(rbind, lapply(out, function(x) cbind(xid = basename(x$file[1]), x[[1]], x[3:6])))
xids = gsub('.RData', '', tstrsplit(merged$xid, '-')[[4]])
reps = as.numeric(tstrsplit(xids, '_')[[2]])
merged$env = factor(reps>8, labels=c('N', 'S'))
merged$rep = reps

ntracks = subset(merged, env=="N")
xgb.load('/Users/noble/Documents/github/cemee/multiwormtracker/male/tp_ngm_xgb_preds')

ntracks = subset(merged, env=="S")
xgb.load('/Users/noble/Documents/github/cemee/multiwormtracker/male/tp_nacl_xgb_preds')

for(i in c('area.F', 'width.F')){
  ntracks[,i] = log(ntracks[,i])
  names(ntracks)[names(ntracks)==i] = paste0('ln.', names(ntracks)[names(ntracks)==i])
}

xpreds = cbind(ntracks[,c('worm', 'xid', 'env', 'rep')], pred = predict(xmod, as.matrix(ntracks[,traits])))
# ggplot(xpreds, aes(pred)) + geom_histogram()
# ggplot(xpreds, aes(pred)) + geom_histogram() + scale_y_log10() + geom_vline(aes(xintercept=0.65))
# xpreds$psex = factor(xpreds$pred>0.65, labels = c('male', 'herm'))
xpreds$psex <- NA
xpreds$psex[xpreds$pred<0.5] = 'male'
xpreds$psex[xpreds$pred>=0.5] = 'herm'
table(xpreds$psex, useNA = 'always')
# remember this is the track frequency across the whole movie - massively overestimates male freq
mf = table(xpreds$rep, xpreds$psex)
mf = data.frame(cbind(mf, mfreq = mf[,2]/(mf[,1]+mf[,2])))
plot(mf$mfreq)

# now annotate tracks in Parsed data and sample for male frequency.
save(xpreds, file = '~/Downloads/tp_G5_ngm.preds.RData')
save(xpreds, file = '~/Downloads/tp_G5_nacl.preds.RData')
