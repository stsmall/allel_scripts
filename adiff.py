#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:38:22 2017
Pairwise Fst, F2, F3, F3, ABBA-BABA, shared doubletons, admixture plots
requires PLINK, ADMIXTURE, vcftools, allel
@author: scott
"""

import allel
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')


def shared_doubletons(ac_subpops):
    """
    """
    # count doubletons
    d2 = {}
    dshare = {}
    for pop in ac_subpops.keys():
        d2[pop] = ac_subpops[pop].is_doubleton(allele=1)
    for x, y in combinations(ac_subpops.keys(), 2):
        s = np.sum([d2[x], d2[y]], axis=0)
        sd = np.nonzero(s == 2)[0].shape[0]
        dshare["{}-{}".format(x, y)] = sd
    print("{}".format(sd))
    return(dshare)


def jsfs_plot(jsfs, pop1, pop2, fold, save):
    """
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    if fold:
        allel.stats.plot_joint_sfs_folded(jsfs, ax=ax)
        title = "jSFS Folded"
    else:
        allel.stats.plot_joint_sfs(jsfs, ax=ax)
        title = "jSFS Polarized"
    fig.suptitle(title, y=1.02)
    ax.set_ylabel('Alternate allele count, {}'.format(pop1))
    ax.set_xlabel('Alternate allele count, {}'.format(pop2))
    if save:
        fig.savefig("jSFS.{}-{}.pdf".format(pop1, pop2), bbox_inches='tight')
    return(None)


def jsfs(ac_subpops, fold=True, save=False):
    """
    """
    jsfsdict = {}
    for x, y in combinations(ac_subpops.keys(), 2):
        ac1 = ac_subpops[x]
        ac2 = ac_subpops[y]
        if fold:
            jsfs = allel.stats.joint_sfs_folded(ac1[:, :2], ac2[:, :2])
        else:
            jsfs = allel.stats.joint_sfs(ac1[:, 1], ac2[:, 1])  # here 1 is alt count
        jsfs_plot(jsfs, x, y, fold, save)
        jsfsdict["{}-{}".format(x, y)] = jsfs
    return(jsfsdict)
