#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:29:50 2017

@author: scott
"""

import allel
import numpy as np
from itertools import combinations
import matplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')


def wcfst(popdict, ac_subpops, var, blenw=10000, downsample=False):
    """Calculates Weir and Cockerham FST
    """
    fstdict = {}
    posdict = {}
    gtdict = {}
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = var.pos[flt]
        geno = var.gt.compress(flt, axis=0)
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
        fst_wc, se_wc, vb_wc, _ =\
            allel.stats.blockwise_weir_cockerham_fst(geno,
                                                     subpops=[pop1_ix,
                                                              pop2_ix],
                                                     blen=blenw, max_allele=1)
        fst_windowed = plot_fst(vb_wc, posflt, blenw, var.chrm)
        fstdict["{}-{}".format(x, y)] = (snp_fst, fst_wc, se_wc, fst_windowed)
        gtdict["{}-{}".format(x, y)] = geno
        posdict["{}-{}".format(x, y)] = posflt
    return(fstdict, gtdict, posdict)


def hdfst(ac_subpops, var, blenw=10000):
    """ Hudson FST
    """
    fstdict = {}
    acdict = {}
    posdict = {}
    pos = var.pos
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt, axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt, axis=0)[:, :2])
        num, dem = allel.stats.hudson_fst(ac1, ac2)
        snp_fst = num / dem
        fst_hd, se_hd, vb_hd, _ =\
            allel.stats.blockwise_hudson_fst(ac1, ac2, blen=blenw)
        fst_windowed = plot_fst(vb_hd, posflt, blenw, var.chrm)
        fstdict["{}-{}".format(x, y)] = (snp_fst, fst_hd, se_hd, fst_windowed)
        posdict["{}-{}".format(x, y)] = posflt
        acdict["{}-{}".format(x, y)] = (ac1, ac2)
    return(fstdict, acdict, posdict)


def plot_fst(vb, pos, blen, chrom):
    """Genome plot of FST along the chrom
    """
    y = vb
    # use the block centres as the X coordinate
    x = allel.stats.moving_statistic(pos,
                                     statistic=lambda v: (v[0] + v[-1])/2,
                                     size=blen)
    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())
    return(x, y)


def pairFST_fx(c, ac_subpops, var, popdict, wc=False):
    """
    """
    if wc:
        fstdict, gtdict, posdict = wcfst(popdict, var)
    else:
        fstdict, acdict, posdict = hdfst(ac_subpops, var)

    return(fstdict)
