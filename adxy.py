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
import matplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')


def plot_dxy(dxy_blocks, blen, pos, chrom):
    """Genome plot of Dxy along the chrom
    """
    y = dxy_blocks
    # use the block centres as the X coordinate
    x = allel.stats.moving_statistic(pos,
                                     statistic=lambda v: (v[0] + v[-1])/2,
                                     size=blen)
    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$D_{XY}$')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())
    return(x, y)


def pairDXY_fx(c, ac_subpops, var, popdict, blenw=10000):
    """Calculates DXY
    """
    dxydict = {}
    pos = var.pos
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt, axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt, axis=0)[:, :2])
        dxy = allel.stats.diversity.mean_pairwise_difference_between(ac1, ac2)
        dxy_blocks = allel.stats.moving_statistic(dxy, statistic=np.mean,
                                                  size=blenw)
        dxy_m, dxy_se, _, _ = jackknife(dxy_blocks)
        dxy_windowed = plot_dxy(dxy_blocks, blenw, posflt, var.chrm)
        dxy["{}-{}".format(x, y)] = (dxy, dxy_m, dxy_se, dxy_windowed)

    return(dxydict)
