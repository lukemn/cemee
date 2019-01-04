#!/usr/bin/env python

import sys, signal
from sys import stderr, stdout, stdin
import os
from os import path
from os.path import basename, dirname
import time
import h5py
import scipy as sp
import pandas as pd
import numpy as np
from itertools import combinations as combn, combinations_with_replacement as combr, product
from limix.stats import indep_pairwise, linear_kinship
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from tqdm import trange, tnrange

# 2-locus epistasis simulations
# simulate interactions across a heritability range (nSim times), 
# test by LR against additive model.
# To declare genome-wide significance, the null is sampled
# at each locus by resampling responses from the
# additive linear model (nBootPerTest times). Empirical 
# p-values can then be calculated from pooled null p-values 
# (at FWER 1/nBootPerTest). 
#
# Output per call:
# (1) nSim x chromosomes, sites, -log10 pvals, h2
# (2) nSim x nBootPerTest null -log10 pvals
#
# NB: at each sampled pair of sites, filtering is applied to ensure the
# presence of all four genotype classes (see tfilter() below)
# and a threshold on missing data is applied.

h2aa    = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25]
maxr2   = 0.5              # LD threshold
Xfile   = sys.argv[1]      # genotype matrix: integer coded, N x P. header assumed
SNPfile = sys.argv[2]      # marker positions: chrom, pos. header assumed
chrom1  = sys.argv[3]      # subset to chrom1/2
chrom2  = sys.argv[4]
ix      = int(sys.argv[5]) # if 0: run initial LD prune & h5 dump. if >0 simulate and prefix output
outpref = 'sim2d_{}{}_{}'.format(chrom1, chrom2, ix)
nSim    = 1000             # per h2 value
nBootPerTest = 100         # parametric bootstraps against additive model per test


def loadData(XF, snpF, chrom=None, startPos=0, endPos=None, r2=0.99, Kr2=0.99):
    '''
    Load genotypes, SNP positions, optionally subsetting on chromosome, position.
    On first load, genome is pruned of local LD (r^2 < Kr2) and an additive GSM (K.txt) is dumped, 
    then genotypes are further pruned (r^2 < r2) and dumped in .h5 format. 
    '''

    def prune(X, snps, r2, win=100, step=50, iters=2):
        '''2-pass local pairwise LD filter on NxP array'''
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


def simY(M, h2fg=[0.1,0,0], h2bg=[0,0,0], Nfg=[1,0,0], Nbg=[0,0,0], dEff='normal'):
    """
    Simulate trait varying additive and epistatic (up to third order) heritability 
    and polygenicity in independent fore- and background components. 
    Pass number of markers to randomly sample (no replacement) in each component 
    and heritability set.

    Effects are drawn from standard normal or exponential dist, with normal residuals.
    Assumes haploid/homozygous diploid genotypes.
    Returns phenotype vector, and dataframe listing causal markers in each set.

    In the special case of 2-locus epistasis, genotypes are filtered on class freq.
    Otherwise, genotypes are continuously sampled until a valid phenotype vector
    can be constructed, which is inefficient for small numbers of causal markers 
    (esp. 3rd order epistasis with only 3 markers).
    """
    
    def scale(X):
        X = X.astype(float)
        X -= X.mean(axis=0)
        X /= X.std(axis=0)
        return(X)
    
    N, P = M.shape
    Vres = 1.0-(sum(h2fg)+sum(h2bg))

    if(P<1 or P<(sum(Nfg)+sum(Nbg))):
        stderr.write('check M dimensions, P={}, Nfg={}, Nbg={}\n'.format(P, Nfg, Nbg))
        exit(1)
    if(N<1):
        stderr.write('check M dimensions, N={}\n'.format(N))
        exit(1)
    if(len(h2fg)<3 or len(h2bg)<3):
        stderr.write('h2 fg/bg must be of length 3\n')
        exit(1)
    if(Vres<0):
        stderr.write('heritability>1\n')
        exit(1)        
    if(Vres==1):
        stderr.write('heritability==0\n')
        exit(1)
    for i in range(3):
        if h2fg[i] and Nfg[i]<(i+1):  
            stderr.write('Insufficient foreground markers {} for heritability {}\n'.format(Nfg, h2fg))
            exit(1)            
        if h2bg[i] and Nbg[i]<(i+1):  
            stderr.write('Insufficient background markers {} for heritability {}\n'.format(Nbg, h2bg))
            exit(1)                       
    
    def addeffs(heritabilities, markers, dEff, y=None, Ix=None): 
        """
        Scale the phenotype vector for a given list of (3) heritabilities 
        and corresponding list of the number of markers for each.
        Returns phenotype vector, Ix vector (trimmed of sampled positions) 
        and list of index/effect arrays for each set.
        If Ix (sampling vector) is passed, then subsequent markers
        are sampled from the beginning of this vector, else all P markers
        are randomly sampled.
        """
        
        def tfilter(array, minN=3, maxNA=5):
            '''
            check two-locus for genotype class frequency, NaN data
            pass = true
            '''

            xtab = pd.crosstab(array[:,0], array[:,1])
            ok=0
            if xtab.shape==(2,2):
                if min(xtab.min()) > minN:
                    if sum(np.isnan(array.ravel())) < maxNA:
                        ok=1
            return(ok)

    
        def effects(dEff, nMarkers):
            if dEff == 'normal': 
                W = sp.randn(nMarkers)
            else:
                W = np.random.exponential(size=nMarkers)
            return(W)
        
        if y is None:
            y = np.zeros(N)
        if Ix is None:
            Ix = sp.random.permutation(P)
        causalMarkers = {}
        
        # additive
        if heritabilities[0]>0:
            i = markers[0]
            I_A = Ix[:i]
            Ix = Ix[i:]
            W_A = effects(dEff, markers[0])
            y_A = sp.dot(M[:,I_A],W_A)
            y += y_A*sp.sqrt(heritabilities[0] / sp.var(y_A))
            causalMarkers['A'] = [I_A, W_A]
    
        # 2nd order epistatic 
        def makeAA(i, Ix):
            '''
            sample and add i effects to genotypes sampled from Ix.
            Return:
              bool flag (true if y is non null)
              y
              index of sites
              effects
            '''

            I_AA = Ix[:i]
            MAA = M[:,I_AA]
            if(i==2):
                while 1:
                    MAA = M[:,I_AA]                
                    if tfilter(MAA):
                        break
                    else:
                        Ix = Ix[1:]
                        I_AA = Ix[:i]
            Ix = Ix[i:]
            k = sp.misc.comb(markers[1], 2, exact=True)
            XAA = np.zeros((N, k))
            for j in range(N):
                gaa = np.outer(MAA[j], MAA[j])
                jj =  gaa[np.triu_indices(gaa.shape[1],1)]
                XAA[j] = jj
            W_AA = effects(dEff, k)
            y_AA = sp.dot(XAA, W_AA)
            if(sp.var(y_AA)):
                return(1, y_AA, I_AA, W_AA)
            else:
                return(0,0,0,0)
            
        if heritabilities[1]>0:
            i = markers[1]
            while 1:
                F, y_AA, I_AA, W_AA = makeAA(i, Ix)
                if F:
                    break
                else:
                    Ix = Ix[1:]  
            Ix = Ix[i:]
            y += y_AA*sp.sqrt(heritabilities[1] / sp.var(y_AA))
            causalMarkers['AA'] = [I_AA, W_AA]
              

        # 3rd order
        def makeAAA(i, Ix):
            I_AAA = Ix[:i]
            MAAA = M[:,I_AAA]
            Ix = Ix[i:]
            k = sp.misc.comb(markers[2], 3, exact=True)
            XAAA = np.zeros((N, k))
            ix3 = combn(np.arange(i), 3)
            for a,j in enumerate(ix3):
                jj = MAAA[:,j[0]]*MAAA[:,j[1]]*MAAA[:,j[2]]
                XAAA[:,a] = jj
                W_AAA = effects(dEff, k)
            y_AAA = sp.dot(XAAA, W_AAA)
            if(sp.var(y_AAA)):
                return(1, y_AAA, I_AAA, W_AAA)
            else:
                return(0,0,0,0)
            
        if heritabilities[2]>0:
            i = markers[2]
            while 1:
                F, y_AAA, I_AAA, W_AAA = makeAAA(i, Ix)
                if F:
                    break
                else:
                    Ix = Ix[1:]
                    if len(Ix)<i:
                        Ix = sp.random.permutation(P)
                        print('failed to find {} suitable markers for AAA in permutation, trying again (components are now potentially non-independent'.format(i))
            Ix = Ix[i:]
            y += y_AAA*sp.sqrt(heritabilities[2] / sp.var(y_AAA)) 
            causalMarkers['AAA'] = [I_AAA, W_AAA]
            
        return(y, Ix, causalMarkers)

    def makeIxDf(ixdic, comp):
        '''make a dataframe of the causal markers in each component, set'''
        ixdf = pd.DataFrame()
        for i in ixdic.keys():
            dfi = pd.DataFrame({'component' : comp, 'set' : i, 'ix' : ixdic[i][0]})
            dfi.sort_values('ix', inplace=True)
            ixdf = pd.concat([ixdf, dfi], ignore_index=True)
        ixdf.reset_index(inplace=True)
        return(ixdf)
    
    y = None
    Ix = None
    ixfg = None
    ixbg = None
    
    # FOREGROUND COMPONENT
    if(sum(h2fg)>0):
        y, Ix, fg = addeffs(h2fg, Nfg, dEff)
        ixfg = makeIxDf(fg, 'fg')
        
    # BACKGROUND COMPONENT
    if(sum(h2bg)>0):
        y, Ix, bg = addeffs(h2bg, Nbg, dEff, y, Ix)
        ixbg = makeIxDf(bg, 'bg')
        
    # normal residual
    y_res = sp.randn(N)
    y += y_res*sp.sqrt(Vres / sp.var(y_res))
    y = scale(y).reshape(N,1)
    
    return(y, ixfg, ixbg)


def hier2d(df, nboot=0):
    '''
    returns:
    1. p-value for full v additive model by least squares, LRT
    2. list of LRT p-values from parametric bootstrap, fitting full model 
    to resampled null responses
    '''
    
    m0 = ols('y ~ x1 + x2', df).fit()
    m1 = ols('y ~ x1 * x2', df).fit()
    t1 = m1.compare_lr_test(m0)[1]
    boots = []    
    Yresp = m0.predict()
    Yresp -= np.mean(Yresp)
    Yresp /= np.std(Yresp)
    N = Yresp.shape[0]    
    dfi = df.copy()
    for i in np.arange(nboot):  
        dfi['y'] = np.random.choice(Yresp, N, replace=1)
        m1i = ols('y ~ x1 * x2', dfi).fit()        
        boots.append(m1i.compare_lr_test(m0)[1])
    return(t1, boots)

def sim2d(X, snp, h2aa, nsim, nboot, outpref):
    
    '''
    generate Y and sample a pair of interacting loci,
    get LRT p-value for full (additive and interaction genetic effects) v 
    additive model, take nperm parametric bootstraps from the additive model
    for pooling to determine genome-wide significance (at FWER 1/nboot).
    '''
    
    out = pd.DataFrame()
    nulls = []
    for i in trange(nsim, desc = 'simulation'):    
        for haa in h2aa:
            y, fg, bg = simY(X, h2fg = [0, haa, 0], Nfg = [0, 2, 0], h2bg = [0, 0, 0], Nbg = [0, 0, 0])
            df = pd.DataFrame(X[:,fg['ix'][fg['set']=='AA']])
            df.columns = map(lambda x: 'x' + str(x+1), df.columns)
            df['y'] = y.ravel()
            stat, boots = hier2d(df, nboot=nboot)
            o = pd.DataFrame({'chrom1' : [snp['chrom'][fg['ix'][0]]],
                             'chrom2' : [snp['chrom'][fg['ix'][1]]],
                             'pos1' : [snp['pos'][fg['ix'][0]]],
                             'pos2' : [snp['pos'][fg['ix'][1]]],
                             'h2' : [haa],
                             'p' : [-sp.log10(stat)],
                             })
            out = pd.concat([out, o])
            nulls += [-sp.log10(boots)]

    out.to_csv(outpref+'_results.txt', index=0, sep='\t')
    nulls = np.array(nulls)
    np.savetxt(outpref+'_nulls.txt.gz', nulls, delimiter='\t', fmt='%.5f')
            
def main():
    X1, snp1, K1 = loadData(Xfile, SNPfile, chrom1, r2=maxr2)
    X2, snp2, K2 = loadData(Xfile, SNPfile, chrom2, r2=maxr2)
    X = np.concatenate([X1, X2], axis=1)
    snp = pd.concat([snp1, snp2], ignore_index=True)[['chrom', 'pos']]
    sim2d(X, snp, h2aa, nSim, nBootPerTest, outpref)

main()

