#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 16:09:58 2017

@author: scott
"""
import allel
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_style('white')
sns.set_style('ticks')


def sfs_plot(c, ac_subpops, save=True, fold=True, scale=True):
    """
    note: should filter on segregating if only using subset of pops
    note: only biallelic if >1 allele
    is_biallelic_01 = ac_seg['all'].is_biallelic_01()[:]
    ac1 = ac_seg['BFM'].compress(is_biallelic_01, axis=0)[:, :2]
    ac2 = ac_seg['AOM'].compress(is_biallelic_01, axis=0)[:, :2]
    """
    sfsdict = {}
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.despine(ax=ax, offset=10)
    for pop in ac_subpops.keys():
        acu = ac_subpops[pop]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('SFS : retaining', np.count_nonzero(flt), 'SNPs')
        # ac1 = allel.AlleleCountsArray(ac_subpops[pop].compress(flt, axis=0)[:, :2])
        ac1 = allel.AlleleCountsArray(ac_subpops[pop].compress(flt, axis=0))
        if fold and scale:
            sfs = allel.sfs_folded_scaled(ac1)
        elif fold and not scale:
            sfs = allel.sfs_folded(ac1)
        elif not fold and not scale:
            sfs = allel.sfs(ac1[:, 1])
        elif not fold and scale:
            sfs = allel.sfs_scaled(ac1[:, 1])
        sfsdict[pop] = sfs
        allel.stats.plot_sfs_folded_scaled(sfsdict[pop], ax=ax, label=pop,
                                           n=ac1.sum(axis=1).max())
    ax.legend()
    ax.set_title('{} Scaled folded site frequency spectra'.format(c))
    ax.set_xlabel('minor allele frequency')
    if save:
        fig.savefig("ScaledSFS-{}.pdf".format(c), bbox_inches='tight')
    return(sfsdict)
