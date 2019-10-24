#!/usr/bin/env python

import sys, signal
from sys import stderr, stdout, stdin
import os
from os import path
from glob import glob
from os.path import basename, dirname
import subprocess
import time
import h5py
import scipy as sp
import pandas as pd
import numpy as np
from limix.stats import Chi2Mixture
from limix.qc import indep_pairwise
from statsmodels.formula.api import ols
from tqdm import trange, tnrange

# 2-locus epistasis simulations
# simulate interactions across a heritability range (nSim times), 
# test by LR against additive model.
# To declare genome-wide significance, the null is sampled
# at each locus by resampling responses from the
# additive linear model (nBootPerTest times) and saving the LR.
# Empirical p-values can then be calculated from pooled null LRs
# when ix=0 (at FWER 1/nBootPerTest). 
#
# Output per call:
# (1) nSim x chromosomes, sites, -log10 pvals, h2
# (2) nSim x nBootPerTest null -log10 pvals
#
# NB: at each sampled pair of sites, filtering is applied to ensure the
# presence of all four genotype classes (see tfilter() below)
# and a threshold on missing data is applied.

h2aa    = [0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2]
maxr2   = 0.5              # LD threshold
Xfile   = sys.argv[1]      # genotype matrix: integer coded, N x P. header assumed
SNPfile = sys.argv[2]      # marker positions: chrom, posi, maf. header assumed
chrom1  = sys.argv[3]      # subset to chrom1/2
chrom2  = sys.argv[4]
ix      = int(sys.argv[5]) # if 0: run initial LD prune & h5 dump. if >0 simulate and prefix output
outpref = 'sim2d_{}{}_{}'.format(chrom1, chrom2, ix)
nSim    = 500              # per h2 value
nBootPerTest = 100         # parametric bootstraps against additive model per test


def loadData(XF, snpF, chrom=None, startPos=0, endPos=None, r2=0.99):
    '''
    Load genotypes, SNP positions, optionally subsetting on chromosome, position.
    On first load, genotypes are pruned of local LD and dumped in .h5 format. 
    '''

    def prune(X, snps, r2, win=100, step=50, iters=2):
        '''local pairwise LD filter on NxP array'''
        for i in range(iters):
            ix = indep_pairwise(X, win, step, r2)
            X = X[:,ix]
            snps = snps.loc[ix]
        return(X, snps)

    if chrom is not None:
        h = '{}.{}.pruned_{}'.format(XF, chrom, r2).replace('.txt', '')
    else:    
        h = XF+'.pruned_{}'.format(r2).replace('.txt', '')
    h5f = h+'.h5'
    snpf = h+'.snps'

    if path.exists(h5f):
        print('loading {}'.format(h))
        snps = pd.read_csv(snpf, sep='\t')
        with h5py.File(h5f, 'r') as hf:
            X = hf["X"][:]            	
    else:
        print('dumping {}'.format(h))
        X = np.loadtxt(XF, delimiter='\t', skiprows=1).T
        snps = pd.read_csv(snpF, sep='\t')
        snps['chrom'] = list(map(str, snps['chrom'].values))
        if chrom is not None:
            q = snps['chrom']==chrom	
            if endPos is not None:
                q = (snps['chrom']==chrom) & (snps['pos']>startPos) & (snps['pos']<endPos) 
            X = X[:,q]
            snps = snps.loc[q]
        m_in = snps.shape[0]
        print('{} (of {}) markers to prune at r^2 < {}'.format(X.shape[1], m_in, r2))
        X, snps = prune(X, snps, r2)
        print('...{} pass'.format(X.shape[1]))
        snps.reset_index(inplace=True)            
        snps.to_csv(snpf, sep='\t')                
        with h5py.File(h5f, 'w') as hf:
            hf.create_dataset("X",  data=X)        

    return(X, snps, h)

def hier2d(df, nboot=0):
    '''
    returns:
    1. p-value for full v additive model by least squares, LRT
    2. adjusted r2 for full model
    3. list of nboot LRT p-values from parametric bootstrap, fitting full model 
    to resampled null responses
    '''
    
    m0 = ols('y ~ x1 + x2', df).fit()
    m1 = ols('y ~ x1 * x2', df).fit()
    r2 = m1.rsquared
    LR1 = m1.compare_lr_test(m0)[0]
    boots = []    
    Yresp = m0.predict()
    Yresp -= np.mean(Yresp)
    Yresp /= np.std(Yresp)
    N = Yresp.shape[0]    
    dfi = df.copy()
    for i in np.arange(nboot):  
        dfi['y'] = np.random.choice(Yresp, N, replace=1)
        m0i = ols('y ~ x1 + x2', dfi).fit()        
        m1i = ols('y ~ x1 * x2', dfi).fit()        
        boots.append(m1i.compare_lr_test(m0i)[0])
    return(LR1, r2, boots)

def sim2d(X, snp, h2aa, nsim, nboot, outpref, h1, h2):
    
    '''
    generate Y (from simAA.R) and sample a pair of interacting loci,
    get LR for full (additive and interaction genetic effects) vs. 
    additive model, take nboot parametric bootstraps from the additive model
    for pooling to determine chisq DoF mixture and p-values (at FWER 1/nboot).
    '''
    
    out = pd.DataFrame()
    nulls = pd.DataFrame()
    for haa in h2aa:
        phef = 'sim_{}_0_{}_2.h5'.format(ix, haa)
        print('Generating phenotypes {}'.format(phef))
        subprocess.call(map(str, ['simAA.R', h1+'.h5', h2+'.h5', 0, haa, 2, nsim, ix]))
        hf = h5py.File(phef, 'r')
        phe = hf['phe']
        markers = pd.DataFrame(hf['ix'][:])
        pk = list(phe.keys())
        for i in trange(nsim, desc = 'simulation {}'.format(haa)):
            df = pd.DataFrame(phe[pk[i]][:])
            LR1, r2, boots = hier2d(df, nboot=nboot)
            j, k = markers.loc()[i]-1
            o = pd.DataFrame({'chrom1' : snp['chrom'][j],
                             'chrom2' : snp['chrom'][k],
                             'pos1' : snp['pos'][j],
                             'pos2' : snp['pos'][k],
                             'maf1' : snp['maf'][j],
                             'maf2' : snp['maf'][k],
                             'h2' : [haa],
                             'r2' : [r2],
                             'LR' : [LR1],
                             })
            out = pd.concat([out, o])
            o = pd.DataFrame({'r2' : [r2]*nboot, 'LR' : boots, 'mafp' : [float(o['maf1'])*float(o['maf2'])]*nboot})
            nulls = pd.concat([nulls, o])

    out.to_csv(outpref+'_results.txt', index=0, sep='\t')
    nulls.to_csv(outpref+'_nulls.txt.gz', index=0, sep='\t', compression='gzip')

def main():
    if ix==0:
        nullf = glob('sim2d*_nulls.txt.gz')
        resf = glob('sim2d*_results.txt')
        print('estimating p-values from {} null files'.format(len(nullf)))        
        pv = getPvals(nullf, resf) 
    else:
        X1, snp1, handle1 = loadData(Xfile, SNPfile, chrom1, r2=maxr2)
        X2, snp2, handle2 = loadData(Xfile, SNPfile, chrom2, r2=maxr2)
        X = np.concatenate([X1, X2], axis=1)
        # indexes from simAA.R are wrt the catted pruned genotype matrices for chrom1,2
        snp = pd.concat([snp1, snp2], ignore_index=True)[['chrom', 'pos', 'maf']]
        sim2d(X, snp, h2aa, nSim, nBootPerTest, outpref, handle1, handle2)

main()

