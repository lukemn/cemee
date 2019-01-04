#!/usr/bin/env python
# coding: utf-8

# In[75]:

import os
from os import path
from os.path import dirname, basename
import sys
import fileinput
import argparse
import h5py
import pandas as pd
import numpy as np
import scipy as sp
import scipy.linalg as splin 
import gc
import gzip
from sklearn import covariance, preprocessing
from limix.stats import indep_pairwise, linear_kinship


# effective number of tests from eigenvalues
# of the ledoit-wolf shrinkage estimator of genotype
# covariances.
# Meff is then simply the number of eigenvalues required to explain
# varex proportion of the variance

maxr2 = 0.9
varex = 0.99
Xfile = "/scratch/lmn3/cemee/sims/2d/gt.txt"
SNPfile = "/scratch/lmn3/cemee/sims/2d/snps"
Xfile ="/scratch/cgsb/rockman/lmn3/G140RILs/plink/v2/CeMEE_v2.0.99.0.05.geno.txt"
SNPfile ="/scratch/cgsb/rockman/lmn3/G140RILs/plink/v2/CeMEE_v2.0.99.0.05.snps"


# In[43]:

def lw_shrink(genotypes):
    '''
    Function to obtain smoothed estimate of the genotype correlation matrix.
    Uses the method proposed by Ledoit and Wolf to estimate the shrinkage parameter alpha.
    Input: genotype matrix.
    Output: smoother correlation matrix and the estimated shrinkage parameter.
    '''
    lw = covariance.LedoitWolf()
    m, n = np.shape(genotypes)
    try:
        fitted = lw.fit(genotypes.T)
        alpha = fitted.shrinkage_
        shrunk_cov = fitted.covariance_
        shrunk_precision = np.mat(np.diag(np.diag(shrunk_cov)**(-.5)))
        shrunk_cor = shrunk_precision * shrunk_cov * shrunk_precision
    except: #Exception for handling case where SNPs in the window are all in perfect LD
        row = np.repeat(1, m)
        shrunk_cor = []
        for i in range(0,m):
            shrunk_cor.append(row)
        shrunk_cor = np.mat(shrunk_cor)
        alpha = 'NA'
    return shrunk_cor, alpha

def find_num_eigs(cov, var_thresh):
    '''
    Function to find the number of eigenvalues required to reach a certain threshold of variance explained.
    '''
    p,p = cov.shape
    eigenv = splin.eigvalsh(cov)
    eigenv[eigenv<0] = 0
    eigenv = np.sort(eigenv).tolist()
    eigenv.reverse()
    running_sum = 0
    counter = 0
    while running_sum < p * var_thresh:
        running_sum += eigenv[counter]
        counter += 1
    return counter


# In[78]:

def loadData(XF, snpF, chrom=None, startPos=0, endPos=None, r2=0.99, Kr2=0.99):
    '''
    Load genotypes, SNP positions, optionally subsetting on chromosome, position.
    On first load, genome is pruned of local LD (r^2 < Kr2) and an additive GSM (K.txt) is dumped, 
    then genotypes are further pruned (r^2 < r2) and dumped to .h5 format. 
    '''

    def prune(X, snps, r2, win=100, step=50, iters=2):
        '''2-pass local pairwise LD filter on NxP array'''
        for i in range(iters):
            ix = indep_pairwise(X, win, step, r2)
            X = X[:,ix]
            snps = snps.loc[ix]
        return(X, snps)

    if chrom is not None:
        chrom = str(chrom)
        h = '{}.{}.pruned_{}'.format(XF, chrom, r2).replace('.txt', '')
    else:
        h = XF+'.pruned_{}'.format(r2).replace('.txt', '')
    h5f = h+'.h5'
    snpf = h+'.snps'

    if path.exists(h5f):
        print('loading {}'.format(h))
        snps = pd.read_csv(snpf, sep='\t')
        K = np.loadtxt(path.join(dirname(XF), 'K.txt'))
        with h5py.File(h5f, 'r') as hf:
            X = hf["X"][:]
    else:
        print('dumping {}'.format(h))
        X = np.loadtxt(XF, delimiter='\t', skiprows=1).T
        snps = pd.read_csv(snpF, sep='\t')
        snps['chrom'] = map(str, snps['chrom'].values)
        X, snps = prune(X, snps, Kr2)
        if not path.exists('K.txt'):
            print('Dumping GSM based on {} (of {}) markers at r^2 < {}'.format(X.shape[1], snps.shape[0], Kr2))
            K = linear_kinship(X, verbose=0)
            np.savetxt(path.join(dirname(XF), 'K.txt'), K)
        else:
            K = np.loadtxt(path.join(dirname(XF), 'K.txt'))
        if chrom is not None:
            q = snps['chrom']==chrom
            if endPos is not None:
                q = (snps['chrom']==chrom) & (snps['pos']>startPos) & (snps['pos']<endPos)
            X = X[:,q]
            snps = snps.loc[q]
        m_in = snps.shape[0]
        X, snps = prune(X, snps, r2)
        print('{} (of {}) markers remain for testing at r^2 < {}'.format(X.shape[1], m_in, r2))
        snps.reset_index(inplace=True)
        snps.to_csv(snpf, sep='\t')
        with h5py.File(h5f, 'w') as hf:
            hf.create_dataset("X",  data=X)
    return(X, snps, K)


# In[79]:

a,b,c=0,0,0
for i in range(6):
    X, snp, K = loadData(Xfile, SNPfile, i+1, r2=maxr2)
    S = preprocessing.scale(X)
    gcov, alpha = lw_shrink(S.T)
    meff_lw = find_num_eigs(gcov, varex)
    a+=meff_lw
    gcov = covariance.empirical_covariance(S, assume_centered=1)
    meff_ep = find_num_eigs(gcov, varex)
    b+=meff_ep
    c+=gcov.shape[0]
    print(i+1, gcov.shape[0], meff_lw, meff_ep)
print(a,b,c)


# In[ ]:



