#!/usr/bin/env Rscript
require(data.table)

ped = commandArgs(trailingOnly = T)[1]
ped = '~/Documents/cemee/v2_manuscript/supplement/CeMEE_v2.MAF0.05.R0.5.ped'
m <- fread(cmd = sprintf('cut -f7- -d" " %s', ped))
m = m[,seq(1, ncol(m), 2),with=F]
fwrite(m, file = file.path(dirname(ped), gsub('.ped', '.txt', basename(ped))), sep = '\t', col.names = F, quote=F)
