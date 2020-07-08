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


def plot_dxy(dxydict, pop2color, pops, chrom, chrsize, save=False):
    """Genome plot of FST along the chrom
    """
    for p in pops:
        f = [key for key in dxydict.keys() if p in key]
        fig, ax = plt.subplots(figsize=(10, 4))
        sns.despine(ax=ax, offset=5)
        title = "{}_Dxy".format(p)
        for pp in f:
            nx = dxydict[pp][2][1]
            x = [(np.sum(i)-1)/2 for i in nx]  # need midpoints
            y = dxydict[pp][2][0]
            ax.plot(x, y, lw=.5, label=pp)
        ax.set_ylabel('$D_{XY}$')
        ax.set_xlabel('Chromosome {} position (bp)'.format(chrom))
        ax.set_xlim(0, chrsize)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        fig.tight_layout()
        if save:
            fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


def pairDxy(c, chrsize, ac_subpops, pos, pop2color, plot=False, blenw=10000, nwindow=100):
    """Calculates DXY
    """
    dxydict = {}
    windlen = int(chrsize / nwindow)
    for x, y in combinations(ac_subpops.keys(), 2):
        # segregating only ?
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print("{} retaining {} SNPs".format("{}-{}".format(x, y),
                                            np.count_nonzero(flt)))
        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                             axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt,
                                                             axis=0)[:, :2])
        # all sites
#        ac1 = ac_subpops[x]
#        ac2 = ac_subpops[y]
#        posflt = pos
        # whole chrom
        dxy = allel.windowed_divergence(posflt, ac1, ac2, size=blenw,
                                        start=1, stop=chrsize)
        dxy_m, dxy_se, *f = jackknife(dxy[0])
        dxy_windowed = allel.windowed_divergence(posflt, ac1, ac2,
                                                 size=windlen, start=1,
                                                 stop=chrsize)
        dxy4plot = (dxy_windowed[0], dxy_windowed[1])
        dxydict["{}-{}".format(x, y)] = (dxy_m, dxy_se, dxy4plot)
    if plot:
        plot_dxy(dxydict, pop2color, list(ac_subpops.keys()), c, chrsize)
    return(dxydict)
