#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:29:50 2017

@author: scott
"""

import allel
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')


def plot_fst(fstdict, pops, chrom, csize, save=False):
    """Genome plot of FST along the chrom
    """
    for p in pops:
        f = [key for key in fstdict.keys() if p in key]
        fig, ax = plt.subplots(figsize=(10, 4))
        sns.despine(ax=ax, offset=5)
        title = "{}_FST".format(p)
        for pp in f:
            ax.plot(fstdict[pp][3][0], fstdict[pp][3][1], lw=.5, label=pp)
        ax.set_ylabel('$F_{ST}$')
        ax.set_xlabel('Chromosome {} position (bp)'.format(chrom))
        ax.set_xlim(0, csize)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        fig.tight_layout()
        if save:
            fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


def fstwindWC(geno, pop1_ix, pop2_ix, pos, blen):
    """
    """
    wc, se, vb, _ =\
        allel.stats.blockwise_weir_cockerham_fst(geno,
                                                 subpops=[pop1_ix, pop2_ix],
                                                 blen=blen, max_allele=1)
    x = allel.stats.moving_statistic(pos,
                                     statistic=lambda v: (v[0] + v[-1])/2,
                                     size=blen)
    return(x, vb)


def fstwindHD(ac1, ac2, pos, blen):
    """
    """
    fst, se, vb, _ = allel.stats.blockwise_hudson_fst(ac1, ac2, blen)
    x = allel.stats.moving_statistic(pos,
                                     statistic=lambda v: (v[0] + v[-1])/2,
                                     size=blen)
    return(x, vb)


def fstchromHD(ac1, ac2, blen):
    """
    """
    fst, se, vb, _ = allel.stats.blockwise_hudson_fst(ac1, ac2, blen)
    return(fst, se)


def fstchromWC(geno, pop1_ix, pop2_ix, blenw):
    """Chromosome average
    """
    wc, se, vb, _ =\
        allel.stats.blockwise_weir_cockerham_fst(geno,
                                                 subpops=[pop1_ix, pop2_ix],
                                                 blen=blenw, max_allele=1)
    return(wc, se)


def wcfst(c, csize, popdict, ac_subpops, var, plot, blenw=10000,
          downsample=False, nwindow=100):
    """Calculates Weir and Cockerham FST
    """
    fstdict = {}
    posdict = {}
    gtdict = {}
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print("{} retaining {} SNPs".format("{}-{}".format(x, y),
                                            np.count_nonzero(flt)))
        posflt = var.pos[flt]
        geno = var.gt.compress(flt, axis=0)
        # TODO: downsample
        if downsample:
            l1 = len(popdict[x])
            l2 = len(popdict[y])
            ds = min([l1, l2])
            if l1 < l2:
                randomsample = np.random.choice(popdict[x], ds, replace=False)
                pop1_ix = randomsample
                pop2_ix = popdict[y]
            elif l2 < l1:
                randomsample = np.random.choice(popdict[y], ds, replace=False)
                pop1_ix = popdict[x]
                pop2_ix = randomsample
            else:
                pass
        else:
            pop1_ix = popdict[x]
            pop2_ix = popdict[y]
        a, b, c = allel.stats.weir_cockerham_fst(geno,
                                                 subpops=[pop1_ix, pop2_ix],
                                                 max_allele=1)
        snp_fst = (a / (a + b + c))[:, 0]  # dist FST
        snp_fst = snp_fst[~np.isnan(snp_fst)]
        # chromosome avg
        fst_wc, se_wc = fstchromWC(geno, pop1_ix, pop2_ix, blenw)
        windlen = int(csize / nwindow)
        fst_windowed = fstwindWC(geno, pop1_ix, pop2_ix, posflt, windlen)
        fstdict["{}-{}".format(x, y)] = (snp_fst, fst_wc, se_wc, fst_windowed)
        gtdict["{}-{}".format(x, y)] = geno
        posdict["{}-{}".format(x, y)] = posflt
    if plot:
        plot_fst(fstdict, list(ac_subpops.keys()), c, csize)
    return(fstdict, gtdict, posdict)


def hdfst(c, csize, ac_subpops, pos, plot, blenw=10000, nwindow=100):
    """ Hudson FST
    """
    fstdict = {}
    acdict = {}
    posdict = {}
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print("{} retaining {} SNPs".format("{}-{}".format(x, y),
                                            np.count_nonzero(flt)))
        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                             axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt,
                                                             axis=0)[:, :2])
        num, dem = allel.stats.hudson_fst(ac1, ac2)
        snp_fst = num / dem
        fst_hd, se_hd = fstchromHD(ac1, ac2, blenw)
        windlen = int(csize / nwindow)
        fst_windowed = fstwindHD(ac1, ac2, posflt, windlen)
        fstdict["{}-{}".format(x, y)] = (snp_fst, fst_hd, se_hd, fst_windowed)
        posdict["{}-{}".format(x, y)] = posflt
        acdict["{}-{}".format(x, y)] = (ac1, ac2)
    if plot:
        plot_fst(fstdict, list(ac_subpops.keys()), c, csize)
    return(fstdict, acdict, posdict)


def pairFST(c, csize, ac_subpops, var, popdict, wc=False, plot=False):
    """
    """
    if wc:
        fstdict, gtdict, posdict = wcfst(c, csize, popdict, ac_subpops, var,
                                         plot)
    else:
        fstdict, acdict, posdict = hdfst(c, csize, ac_subpops, var.pos, plot)

    return(fstdict)
