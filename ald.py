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
import scipy
from collections import defaultdict
sns.set_style('white')
sns.set_style('ticks')


def plot_lddecay(c, lddict, xmax, xbin=100, save=True):
    """
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    chrm = c
    title = "LD decay"
    for pop in lddict.keys():
        diff = lddict[pop][0]
        r2 = lddict[pop][1]
        s, e, _ = scipy.stats.binned_statistic(diff[:], r2[:],
                                               statistic=np.nanmean,
                                               bins=np.arange(0, xmax + xbin,
                                                              xbin))
        y = s
        x = (e[:-1] + e[1:]) / 2
        ax.plot(x, y, lw=1.5, label=pop)
        ax.set_ylabel("r2")
        ax.set_xlabel('Chromosome {} position (bp)'.format(chrm))
        ax.set_xlim(0, xmax)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        fig.tight_layout()
        if save:
            fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


def plot_lddict(lddict, xmax=10000, xbin=100, save=True):
    """
    """
    tmpdict = defaultdict(list)
    for c in lddict.keys():
        for p in lddict[c].keys():
            tmpdict[p].append(lddict[c][0])
    diffdict = {}
    for p in tmpdict.keys():
        x = np.concatenate(tmpdict[p])
        diffdict[p] = x[~np.isnan(x)]
    tmpdict = defaultdict(list)
    for c in lddict.keys():
        for p in lddict[c].keys():
            tmpdict[p].append(lddict[c][1])
    r2dict = {}
    for p in tmpdict.keys():
        x = np.concatenate(tmpdict[p])
        r2dict[p] = x[~np.isnan(x)]
    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    title = "LD decay all chromosomes"
    for pop in lddict.keys():
        diff = diffdict[pop]
        r2 = r2dict[pop]
        s, e, _ = scipy.stats.binned_statistic(diff[:], r2[:],
                                               statistic=np.nanmean,
                                               bins=np.arange(0, xmax + xbin,
                                                              xbin))
        y = s
        x = (e[:-1] + e[1:]) / 2
        ax.plot(x, y, lw=.5, label=pop)
        ax.set_ylabel("r2")
        ax.set_xlabel('Distance position (bp)')
        ax.set_xlim(0, xmax)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        fig.tight_layout()
        if save:
            fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


def pdiff(x):
    """
    """
    x = np.subtract.outer(x, x)
    o, t = np.tril_indices_from(x, k=-1)
    m = sorted(zip(o, t), key=lambda x: x[1])
    return(x[list(zip(*m))])


def ld_decay(c, chrsize, ac_subpops, popdict, pop2color, var, min_maf=0.1,
             xmax=6000):
    """
    """
    lddict = {}
    fig, ax = plt.subplots()
    for x in ac_subpops.keys():
        acu = ac_subpops[x]
        pos = var.pos
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        pos = pos[flt]
        ac = allel.AlleleCountsArray(ac_subpops[x].compress(flt, axis=0)[:, :2])
        gt = var.gt.compress(flt, axis=0)[:, popdict[x]]
        af = ac.to_frequencies()
        flt = (af[:, :2].min(axis=1) > min_maf)
        pos = pos[flt]
        gt = gt.compress(flt, axis=0)
        gn = gt.to_n_alt()
        print("calc r2...")
        r = allel.stats.rogers_huff_r(gn) ** 2
        print("calc pdist...")
        dist = pdiff(pos)
        xmax_diff = dist <= xmax
        r2 = r[xmax_diff]
        diff = dist[xmax_diff]
        lddict[x] = (diff, r2)
    plot_lddecay(c, lddict, xmax)
    return(lddict)
