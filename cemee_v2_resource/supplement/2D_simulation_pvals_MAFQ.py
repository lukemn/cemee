#!/usr/bin/env python

import sys, signal
from sys import stderr, stdout, stdin
import os
from os import path
from glob import glob
from os.path import basename, dirname
import scipy as sp
import pandas as pd
import numpy as np
from limix.stats import Chi2Mixture

# 2-locus epistasis simulations
# estimate empirical p-values from pooled null LRs, 
# stratifed by joint MAF.
# pass mq MAF quantile to process and nq total number of quantiles
mq, nq = list(map(int, sys.argv[1:]))

def getPvals(nullf, resf, mq, nquantile):
    '''
    pool null LLRs and estimate X2 components (0, >0) for each parameter set
    stratify by joint MAF quantiles
    '''
    def LR2p(name, group, res):
        x2m = Chi2Mixture(tol=1e-3)
        x2m.estimate_chi2mixture(np.array(group['LR'].astype(float)))
        pv = x2m.sf(res['LR'].astype(float))
        res = res.assign(p = -sp.log10(pv))
        res.to_csv('sim2d_maf{}.pv.txt'.format(name), sep='\t')

    qlabs = np.arange(nquantile)+1
    nulls = (pd.read_csv(f, sep = '\t', dtype=object) for f in nullf)
    res = (pd.read_csv(f, sep = '\t', dtype=object) for f in resf)
    nulldf = pd.concat(nulls, ignore_index=True)
    resdf = pd.concat(res, ignore_index=True)
    resdf['mafq'] = pd.qcut(pd.to_numeric(resdf['maf1'])*pd.to_numeric(resdf['maf2']), nquantile, labels=qlabs)
    nulldf['mafq'] = pd.qcut(pd.to_numeric(nulldf['mafp']), nquantile, labels=qlabs)
    nullsub = nulldf.loc()[nulldf['mafq']==mq]
    ressub = resdf.loc()[resdf['mafq']==mq]
    print('processing {} null {} alt LRs for MAF quantile {}/{}'.format(nullsub.shape[0], ressub.shape[0], mq, nquantile))
    LR2p(mq, nullsub, ressub)

nullf = glob('sim2d*_nulls.txt.gz')
resf = glob('sim2d*_results.txt')
print('estimating p-values from {} null files'.format(len(nullf)))
getPvals(nullf, resf, mq, nq)


