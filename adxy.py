#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:39:58 2017

@author: scott
"""
from autil import jackknife
import allel
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')


def windDxy(dxy, blen, pos):
    """Genome plot of Dxy along the chrom
    """
    dxy_blocks = allel.stats.moving_statistic(dxy, statistic=np.nanmean,
                                              size=blen)
    y = dxy_blocks
    # use the block centres as the X coordinate
    x = allel.stats.moving_statistic(pos,
                                     statistic=lambda v: (v[0] + v[-1])/2,
                                     size=blen)
    return(x, y)


def plot_dxy(dxydict, pops, chrom, save=False):
    """Genome plot of FST along the chrom
    """
    for p in pops:
        f = [key for key in dxydict.keys() if p in key]
        fig, ax = plt.subplots(figsize=(10, 4))
        sns.despine(ax=ax, offset=5)
        m = []
        title = "{}_FST".format(p)
        for pp in f:
            m.extend(dxydict[pp][3][0])
            ax.plot(dxydict[pp][3][0], dxydict[pp][3][1], lw=.5, label=pp)
        ax.set_ylabel('$F_{ST}$')
        ax.set_xlabel('Chromosome {} position (bp)'.format(chrom))
        ax.set_xlim(0, max(m))
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        fig.tight_layout()
        if save:
            fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


def pairDxy(c, ac_subpops, var, popdict, blenw=10000, plot=False):
    """Calculates DXY
    """
    dxydict = {}
    pos = var.pos
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                             axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt,
                                                             axis=0)[:, :2])
        dxy = allel.stats.diversity.mean_pairwise_difference_between(ac1, ac2)
        dxy_blocks = allel.stats.moving_statistic(dxy, statistic=np.nanmean,
                                                  size=blenw)
        dxy_m, dxy_se, _, _ = jackknife(dxy_blocks)
        dxy_windowed = windDxy(dxy, int(blenw/5), posflt)
        dxydict["{}-{}".format(x, y)] = (dxy, dxy_m, dxy_se, dxy_windowed)
    if plot:
        plot_dxy(dxydict, list(ac_subpops.keys()), var.chrm)
    return(dxydict)
