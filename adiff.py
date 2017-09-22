#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:38:22 2017
Pairwise Fst, F2, F3, F3, ABBA-BABA, shared doubletons, admixture plots
requires PLINK, ADMIXTURE, vcftools, allel
@author: scott
"""

import allel
import seaborn as sns
import numpy as np
import pandas as pd
from autil import jackknife
from asfs import jsfs_fx
sns.set_style('white')
sns.set_style('ticks')


def doubletons_fx(ac_subpopscat, subpops):
    """
    double: dict(list)
        PNG: ([False, False, True ...
    """
    # count doubletons
    d2 = {}
    dshare = {}
    for pop in subpops.keys():
        d2[pop] = ac_subpopscat[pop].count_doubleton(allele=1)
    for pop in d2.keys():
        if "_" in pop:
            if "Africa" in pop:
                dshare[pop] = d2[pop.split("_")[0]] + d2["Mali"]\
                    + d2["Kenya"] - d2[pop]
            else:
                dshare[pop] = d2[pop.split("_")[0]] + d2[pop.split("_")[1]]\
                    - d2[pop]
        else:
            pass
    print("{}".format(d2))
    print("{}".format(dshare))
    return(d2, dshare)


def allelecount_pair(ac_subpopsdict, ac_subpopscat, chrcat_gt, subpops,
                     pairlist, jsfs_r, chrpos):
    """
    """
    # jsfs
    if jsfs_r:
        jsfs = allel.stats.joint_sfs(ac1[:, 1], ac2[:, 1])
        jsfsf = allel.stats.joint_sfs_folded(ac1[:, :2], ac2[:, :2])
        jsfs_fx(jsfs, jsfsf, pop1, pop2)

    return(None)
