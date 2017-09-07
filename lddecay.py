#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:51:31 2017

@author: scott
"""
import allel
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import bcolz
import scipy
sns.set_style('white')
sns.set_style('ticks')


def snp_ascertainment(chrdict, nchr, ac, min_maf, popidx):
    """
    """
    pos = chrdict[nchr].positions
    geno = chrdict[nchr].genotypes
    start = 0
    stop = max(pos)
    # ascertainment
    af = ac.to_frequencies()
    loc_asc = (ac.max_allele() < 2) & (af[:, :2].min(axis=1) > min_maf)
    loc_region = np.zeros(pos.size, dtype='b1')
    loc_region[pos.locate_range(start, stop)] = True
    loc_asc &= loc_region
    gt = geno.subset(loc_asc, popidx)
    gn = gt.to_n_alt()

    return(gn, pos[loc_asc])


def pdiff(x):
    """
    """
    x = np.subtract.outer(x, x)
    o, t = np.tril_indices_from(x, k=-1)
    m = sorted(zip(o, t), key=lambda x: x[1])
    return(x[list(zip(*m))])


def plot_ld_decay(chrdict, ac_subpopsdict, subpops, min_maf, xmax, xbin,
                  chrlist, label=None, ax=None):
    """
    """
    pop_color = {"PNG": 'red',
                 "Haiti": 'blue',
                 "Mali": 'green',
                 "Kenya": 'purple'}
    fig, ax = plt.subplots()
    poplist = ["Haiti", "Mali", "Kenya", "PNG"]
    chromlist = ["Wb_Chr3_1"]
    # chromlist = chrlist
    for nchr in chromlist:
        for pop in poplist:
            popidx = subpops[pop]
            ac = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop][:])
            # get genotypes
            gn, pos = snp_ascertainment(chrdict, nchr, ac, min_maf, popidx)
            r = allel.stats.rogers_huff_r(gn) ** 2
            dist = pdiff(pos)
            xmax_diff = dist <= xmax
            r2 = r[xmax_diff]
            diff = dist[xmax_diff]
            s, e, _ = scipy.stats.binned_statistic(diff[:], r2[:],
                                                   statistic=np.nanmean,
                                                   bins=np.arange(0,
                                                                  xmax + xbin,
                                                                  xbin))
            if ax is None:
                fig, ax = plt.subplots()
            x = (e[:-1] + e[1:]) / 2
            y = s
            ax.plot(x, y, marker=' ', linestyle='-', color=pop_color[pop],
                    lw=2, label=pop)
            ax.legend(fontsize='large')
            ax.set_ylim(0, .7)
            # ax.set_xscale('log')
            sns.set_style('ticks')
            sns.despine(offset=5)
            fig.savefig("LDdecay.{}.pdf".format(nchr), bbox_inches='tight')
    return(None)
