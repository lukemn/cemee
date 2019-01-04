I#!/usr/bin/env python

# Single marker/single trait simulations fitting set-based 2-variance component LMM:
# a focal genetic similarity matrix (R) for each window and a genomic similarity matrix (K). 
# Ideally, these should be independent, but I'm ignoring proximal contamination 
# here so we don't have to recalculate K each step.
# Heritability of a single focal marker and a polygenic component (of equal variance to 
# focal marker, but spread over 100, 500, 1000 markers) is varied.
# see https://www.ncbi.nlm.nih.gov/pubmed/26076425
# https://limix.readthedocs.io/en/stable/mtset.html

# Pass: 
# (1) full genotype matrix: integer coded, markers x lines, header, tsv
# (2) markers positions: [chrom, pos], chroms coded as integers 1-6; header, tsv
# (3) numeric index for output.

# On each call a single chromosome is sampled, then nsim simulations of each 
# scenario are carried out. mtSet returns the log likelihood ratio under H1 for each window,
# then the H0 is sampled by permuting strains in R.

# Outputs (per call, results are appended each simulation; all tsv):
# (1) full table of H1 LLRs/set across simulation grid
# (2) matched null values from permutation of R (default 100 samples/window/simulation)
# (3) the simulated focal sites, and some summary stats (inc. single marker LMM p-values).

# p-values are calculated from null LLRs by estimating the appropriate ratio and 
# degrees of freedom of the 0 : >0 mixed chisq distribution for the LRT
# (see Listgarten, 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673214/)
# If index <0, p-vals are estimated from the parametric null, pooling LLRs 
# from each scenario, and written to mtsetSim.pv.txt. 

import sys, signal
from sys import stderr, stdout, stdin
from os import path
from glob import glob
import scipy as sp
import scipy.linalg
import pandas as pd
import numpy as np
import h5py
from tqdm import tnrange, trange
from limix.qtl import qtl_test_lmm
from limix.stats import linear_kinship, Chi2mixture, indep_pairwise
from limix.util import sets_from_bim
from limix.mtset import MTSet

h2focal         = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15]
Nbgs            = [100, 500, 1000]
chroms          = range(1,7)
maxr2           = 0.98            # LD pruning threshold
window, step    = 25000, 25000    # set window, setp size in bp
minVarPerWindow = 3               # minimum number of variants per set
nsim            = 20              # number of simulations per h2/Nbg scenario

Xfile = sys.argv[1]            # genotypes. PxN
SNPfile = sys.argv[2]          # markers. chrom, pos
ix = int(sys.argv[3])          # 0: initial LD filter & h5 dump. <0: empirical pvals
if ix <0 and len(sys.argv)==5: # optional glob pattern to subset pvals by scenario
    hnull = sys.argv[4]
else:
    hnull = "*"


def loadData(XF, snpF, chrom=None, startPos=0, endPos=None, r2=0.99, Kr2=0.99):
    '''
    Load genotypes, SNP positions, optionally subsetting on chromosome, position.
    On first load, genome is pruned of local LD (r^2 < Kr2) and an additive GSM (K.txt) is dumped for reuse, 
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
        h = '{}.{}.pruned_{}'.format(XF, chrom, r2).replace('.txt', '')
    else:    
        h = XF+'.pruned_{}'.format(r2).replace('.txt', '')
    h5f = h+'.h5'
    snpf = h+'.snps'

    if path.exists(h5f):
        print('loading {}'.format(h))
        snps = pd.read_csv(snpf, sep='\t')
        K = np.loadtxt('K.txt')
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
            np.savetxt('K.txt', K) 
        else:
            K = np.loadtxt('K.txt')
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

def scale(X):
    X = X.astype(float)
    X -= X.mean(axis=0)
    X /= X.std(axis=0)
    return(X)
    
def simY(Vfg, Vbg, M, Nfg=1, Nbg=1, I_fg=None, I_bg=None, dEff='normal'):
    """
    Simulate trait varying additive heritability and polygenicity in two components.
    Pass index array for large-effect (foreground) and small-effect (background) causal markers, 
    or number to randomly sample. Effects are drawn from standard normal or exponential dist.
    Returns phenotype vector, and positions of the causal markers.
    """
    N, P = M.shape
    Vres = 1.0-(Vfg+Vbg)

    Ix = sp.random.permutation(P)        
    if I_fg is not None:
        Nfg = len(I_fg)
    if I_bg is not None:
        Nbg = len(I_bg)
    if I_fg is not None and I_bg is None:
        Ix = Ix[~np.in1d(Ix, I_fg)]
    if(Nfg + Nbg >= P):
        Nbg = P-(Nfg+1)
    
    # causal loci effect distribution (N or exp), N error
    if dEff == 'normal': 
        W_fg  = sp.randn(Nfg)
        W_bg  = sp.randn(Nbg)
    else:
        W_fg  = np.random.exponential(size=Nfg)
        W_bg  = np.random.exponential(size=Nbg)
    y_resid = sp.randn(N)

    # focal snp(s)
    if I_fg is None:
        I_fg = Ix[:Nfg]
    y_fg = sp.dot(M[:,I_fg],W_fg)
    y = y_fg*sp.sqrt(Vfg / sp.var(y_fg))

    # polygenic component
    if Vbg>0:
        if I_bg is None:
            I_bg = Ix[Nfg:Nfg+Nbg]
        y_bg = sp.dot(M[:,I_bg],W_bg)
        y += y_bg*sp.sqrt(Vbg / sp.var(y_bg))
    
    # residual
    y_res = sp.randn(N)
    y += y_res*sp.sqrt(Vres / sp.var(y_res))
    y = scale(y).reshape(N,1)
    
    ixfg = pd.DataFrame({'fg' : I_fg, 'w' : W_fg})
    ixfg.sort_values('fg', inplace=True)
    ixfg.reset_index(inplace=True)
    ixbg = pd.DataFrame({'bg' : I_bg, 'w' : W_bg})
    ixbg.sort_values('bg', inplace=True)
    ixbg.reset_index(inplace=True)    
    
    return(y, ixfg, ixbg)
    
def fitMtset(X, bim, sets, Y, K, ix, traitNames=None, rankR=1, n_perms=10, v=0):
    '''
    fit the mtSet test (default rank=1 for R), 
    and run n_perms permutations of lines within R.
    '''
    def getVariants(set, X):
        q = (bim['chrom']==set['chrom']) & (bim['pos']>=set['start']) & (bim['pos']<=set['end'])
        Xr = X[:,q]
        Xr = Xr[:,~np.isnan(Xr).all(axis=0)]
        Xr = scale(Xr)
        return(Xr)

    mtSet = MTSet(Y=Y, R=K, traitID=traitNames, rank=rankR)
    res_null = mtSet.fitNull()
    n_wnds = sets.shape[0]
    
    errc=0
    LLR = sp.zeros(n_wnds) # test stats
    iter_sets = sets.iterrows()
    for wnd_i in trange(n_wnds, desc='alt', leave = False, disable=not v):            
        _, set_i = next(iter_sets)
        Xr = getVariants(set_i, X)
        try:
            RV = mtSet.optimize(Xr)
            LLR[wnd_i] = RV['LLR'][0]
        except ValueError:
            LLR[wnd_i] = 1
            errc+=1
    if errc: stderr.write('{} models failed\n'.format(errc))

    LLR_null = []
    for perm_i in trange(n_perms, desc='perms', leave = False, disable=not v):
        perm_idxs = sp.random.permutation(Y.shape[0])
        iter_sets = sets.iterrows()
        for wnd_i in range(n_wnds):
            _, set_i = next(iter_sets)
            Xr = getVariants(set_i, X)
            Xr = Xr[perm_idxs, :]
            try:
                RV = mtSet.optimize(Xr)
                LLR_null.append(RV['LLR'][0])
            except ValueError:
                pass
    return(LLR, LLR_null)

def getPvals(nullf, resf):
    '''
    pool null LLRs and estimate X2 components (0, >0) for each parameter set
    write by group (heritability ~ Nbg)
    '''
    def LR2p(name, group, res):
        c2m = Chi2mixture(tol=1e-3)
        c2m.estimate_chi2mixture(np.array(group['LLR'].astype(float)))
        resi = res.get_group(name)
        pv = c2m.sf(resi['LLR'].astype(float))
        resi['p'] = -sp.log10(pv)
        resi.to_csv('mtsetSim_{}_{}.pv.txt'.format(name[0], name[1]), sep='\t')        
     
    nulls = (pd.read_csv(f, sep = '\t', dtype=object) for f in nullf)
    res = (pd.read_csv(f, sep = '\t', dtype=object) for f in resf)    
    nulldf = pd.concat(nulls, ignore_index=True)
    nulldf = nulldf.loc[nulldf['h2']!='h2']
    resdf = pd.concat(res, ignore_index=True)
    resdf = resdf.loc[resdf['h2']!='h2']
    nulls = nulldf.groupby(['h2', 'nbg'])
    res = resdf.groupby(['h2', 'nbg'])
    for name, group in nulls:
        print name, group.shape[0]
        LR2p(name, group, res)
    
def simulateAway(window, step, n_perms=100, v=0, init=False):

    # filter and dump, sum the sets for given LD/window settings and exit
    if init:
        nsets = 0
        for chromi in chroms:
            X, snps, K = loadData(Xfile, SNPfile, chrom = str(chromi), r2=maxr2)
            snps['chrom'] = map(str, snps['chrom'].values)
            sets = sets_from_bim(snps, size=window, step=step, minSnps=minVarPerWindow)
            nsets += len(sets)
        print(nsets)
        exit()

    # pick a chromosome
    chromi = np.random.choice(chroms)
    X, snps, K = loadData(Xfile, SNPfile, chrom = chromi, r2=maxr2)
    snps['chrom'] = map(str, snps['chrom'].values)
    sets = sets_from_bim(snps, size=window, step=step, minSnps=minVarPerWindow)
    
    altf = 'mtset.sim.{}.out.txt'.format(ix)
    nullf = 'mtset.sim.{}.null.txt'.format(ix)
    sitesf = 'mtset.sim.{}.sites.txt'.format(ix)
    writeHeader=True
    if path.exists(altf):
        writeHeader=False

    for i in trange(nsim, desc = 'simulations'):    
        for hi in h2focal:
            for nbg in Nbgs:
                # pick a single marker (allowed outside of a set)
                snpi = [np.random.randint(snps.shape[0])]
                q="chrom == '{}' & start >= {} & end <= {}".format(snps['chrom'].values[snpi][0], snps['pos'].values[snpi][0]-window, snps['pos'].values[snpi][0]+window)
                # generate Y and get positions of all causal markers
                Y, fg, bg = simY(hi, hi, X, Nbg = nbg, I_fg = snpi)
                # fit the model, sample the null
                llr, llr_null = fitMtset(X, snps, sets, Y, K, ix, v=v, n_perms=n_perms)
                        
                sets['fg'] = sets['setid'].isin(sets.query(q)['setid'])
                sets['h2'] = hi
                sets['nbg'] = nbg
                sets['LLR'] = llr                                    
                sets['sim'] = '{}.{}'.format(ix, i)
                with open(altf, 'a') as of:
                    sets.to_csv(of, index=False, header=writeHeader, sep='\t')
                
                null = pd.DataFrame({'LLR' : llr_null})
                null['h2'] = hi
                null['nbg'] = nbg
                with open(nullf, 'a') as of:
                    null.to_csv(of, index=False, header=writeHeader, sep='\t')

                fgsnps = snps.loc[fg['fg']]
                fgsnps['w'] = fg['w'].values
                fgsnps['r2'] = sp.stats.pearsonr(Y,X[:,snpi])[0]**2
                fgsnps['h2'] = hi
                fgsnps['nbg'] = nbg
                fgsnps['lmm_p'] = -sp.log10(qtl_test_lmm(X[:,snpi], Y, K, test='lrt').pvalues.ravel()[0])
                fgsnps['sim'] = '{}.{}'.format(ix, i)
                with open(sitesf, 'a') as of:
                    fgsnps.to_csv(of, index=False, header=writeHeader, sep='\t')


def main():
    if ix < 0:
        # should parallelize with dask
        # for now can be done by awk, and passing hnull prefix
        nullf = glob('{}.null.txt'.format(hnull))
        resf = glob('{}.out.txt'.format(hnull))
        print('estimating p-values from {} null files'.format(len(nullf)))        
        pv = getPvals(nullf, resf)
    elif ix==0:
        simulateAway(window, step, init=True)
    else:
        simulateAway(window, step)

main()

