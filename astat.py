#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:30:34 2017
Assorted functions for using scikit allel
@author: scott
"""
import allel
import numpy as np
import bisect
from allel.util import asarray_ndim, check_dim0_aligned, ensure_dim1_aligned
from itertools import combinations
from aplot import plotmiss


def fixdiff(ac_subpops, popdict):
    """
    """
    diffalleles = {}
    for x, y in combinations(popdict.keys(), 2):
        ac1 = ac_subpops[x]
        ac2 = ac_subpops[y]
        loc_df = allel.locate_fixed_differences(ac1, ac2)
        diffalleles["{}-{}".format(x, y)] = sum(loc_df)
    return(diffalleles)


def privalleles(ac_subpops, popdict):
    """
    """
    acs = ac_subpops.values()
    acs = [asarray_ndim(ac, 2) for ac in acs]
    check_dim0_aligned(*acs)
    acs = ensure_dim1_aligned(*acs)
    pvalleles = np.zeros(len(popdict.keys()), dtype="int32")
    pac = np.dstack(acs)
    for s in range(pac.shape[0]):
        x = np.sum(pac[s], axis=1)
        if np.any(x == 1):
            p = np.where(x == 1)[0]
            p2 = np.where(pac[s][p] == 1)[0]
            pvalleles[p2] += 1
    pva = {}
    for i, pop in enumerate(ac_subpops.keys()):
        pva[pop] = pvalleles[i]
    return(pva)


def misspos(chrm, callset, pc, samples, window_size=10000, title=None,
            saved=False):
    """
    """
#    chrm = chrm.decode("utf-8")
    chrom = callset['variants/CHROM']
    chrom_mask = np.where(chrom[:] == chrm)
    pos = callset['variants/POS']
    p = pos[:][chrom_mask]
    varpos = allel.SortedIndex(p)
    bins = np.arange(0, varpos.max(), window_size)
    # use window midpoints as x coordinate
    x = bins
    miss_site = pc[:][chrom_mask]
    yy = []
    for i, j in enumerate(x):
        try:
            left = bisect.bisect_left(varpos, j)
            right = bisect.bisect_left(varpos, x[i+1]) - 1
            yy.append(np.mean(miss_site[left:right]))
        except:
            yy.append(0)
    y = np.array(yy)
    plotmiss(x, y/samples, title, chrm, saved)
