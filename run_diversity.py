#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 18:15:38 2017

@author: scott
"""

import allel
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_style('white')
sns.set_style('ticks')


def sfs_fx(ac_subpopscat):
    """
    note: should filter on segregating if only using subset of pops
    note: only biallelic if >1 allele
    is_biallelic_01 = ac_seg['all'].is_biallelic_01()[:]
    ac1 = ac_seg['BFM'].compress(is_biallelic_01, axis=0)[:, :2]
    ac2 = ac_seg['AOM'].compress(is_biallelic_01, axis=0)[:, :2]
    """
    ac1 = ac_subpopscat["PNG"]
    ac2 = ac_subpopscat["Haiti"]
    ac3 = ac_subpopscat["Mali"]
    ac4 = ac_subpopscat["Kenya"]
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.despine(ax=ax, offset=10)
    sfs1 = allel.stats.sfs_folded_scaled(ac1)
    allel.stats.plot_sfs_folded_scaled(sfs1, ax=ax, label='PNG',
                                       n=ac1.sum(axis=1).max())
    sfs2 = allel.stats.sfs_folded_scaled(ac2)
    allel.stats.plot_sfs_folded_scaled(sfs2, ax=ax, label='Haiti',
                                       n=ac2.sum(axis=1).max())
    sfs3 = allel.stats.sfs_folded_scaled(ac3)
    allel.stats.plot_sfs_folded_scaled(sfs3, ax=ax, label='Mali',
                                       n=ac3.sum(axis=1).max())
    sfs4 = allel.stats.sfs_folded_scaled(ac4)
    allel.stats.plot_sfs_folded_scaled(sfs4, ax=ax, label='Kenya',
                                       n=ac4.sum(axis=1).max())
    ax.legend()
    ax.set_title('Scaled folded site frequency spectra')
    ax.set_xlabel('minor allele frequency')
    return(None)


def tajd(ac_subpops):
    """
    """
    for pair in pairlist:
        pop1 = pair.split("_")[0]
        pop2 = pair.split("_")[1]
        ac_seg = ac_subpops.compress(genotypes_seg[pair])
    pos = variants_seg['POS'][:]
    windows = allel.stats.moving_statistic(pos,
                                           statistic=lambda v: [v[0],
                                                                v[-1]],
                                                                size=2000)
    x = np.asarray(windows).mean(axis=1)

    # compute Tajima's D
    y1, _, _ = allel.stats.windowed_tajima_d(pos, ac_seg['BFM'][:],
                                             windows=windows)
    y2, _, _ = allel.stats.windowed_tajima_d(pos, ac_seg['AOM'][:],
                                             windows=windows)

    # plot
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y1, lw=.5, label='BFM')
    ax.plot(x, y2, lw=.5, label='AOM')
    ax.set_ylabel("Tajima's $D$")
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    return(None)


def pi():
    """
    """
    return(None)


def ld():
    """
    """
    return(None)
