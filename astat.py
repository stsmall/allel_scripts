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
from collections import defaultdict
from itertools import combinations
import aplot as ap


def fixdiff(ac_subpops, popdict):
    """
    """
    diffalleles = {}
    for x, y in combinations(popdict.keys(), 2):
        ac1 = ac_subpops[x]
        ac2 = ac_subpops[y]
        loc_df = allel.locate_fixed_differences(ac1, ac2)
        diff = sum(loc_df)
        diffalleles["{}-{}".format(x, y)] = diff
        print("{}-{}:{}".format(x, y, diff))
    return(diffalleles)


def privalleles(ac_subpops, popdict):
    """try is_seg on each subpop. compare all subpops and count if only True in
    1 pop
    """
    segs = []
    pops = []
    for pop in ac_subpops.keys():
        ac = ac_subpops[pop].is_segregating()
        segs.append(ac)
        pops.append(pop)
    a = np.vstack(segs)  # where is seg
    b = a.sum(axis=0)  # sum bool
    c = np.where(b == 1)[0]  # find where only 1 is true
    d = a[:, c]  # return array of private
    for p in range(a.shape[0]):
        print("{}:{}".format(pops[p], np.count_nonzero(d[p, :])))
    return(None)


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
    ap.plotmiss(x, y/samples, title, chrm, saved)


def gtstats(calls, pop2color, n_variants):
    """
    """
    gtd = allel.GenotypeDaskArray(calls['calldata/GT'])
    pc_missing = gtd.count_missing(axis=0)[:].compute()  # per sample
    miss = gtd.count_missing(axis=1)[:].compute()
    pc_het = gtd.count_het(axis=0)[:].compute()  # per sample
    dep = calls['calldata/DP']
    dp = np.mean(dep[:, :], axis=0)
    ap.plotstats(pc_het/n_variants, 'Heterozygous', pop2color)
    ap.plotstats(pc_missing/n_variants, 'Missing', pop2color)
    ap.plotstats(dp, 'Depth', pop2color)
    return(miss)


def chrstats(chrlist, calls, miss, nsamples):
    """
    """
    # Chromosome GT Stats
    for c in chrlist:
        ap.plotvars(c, calls)
    for c in chrlist:
        misspos(c, calls, miss, nsamples, window_size=10000, saved=False)


def popstats(chrlist, meta, popdict, var):
    """
    """
    chrstatdict = defaultdict(dict)
    for c in chrlist:
        var.geno(c, meta)
        print("\nStats for Chromosome {}\n".format(c))
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
        # population stats: SNP, singleton, doubleton, #het inds, #homref, alt
        for pop in popdict.keys():
            seg = ac_subpops[pop].count_segregating()
            sing = ac_subpops[pop].count_singleton()
            doub = ac_subpops[pop].count_doubleton()
            print("{} SNPs, singleton, doubleton: {} {} {}".format(pop, seg,
                                                                   sing, doub))
            gt_subpop = var.gt.take(popdict[pop], axis=1)
            het = gt_subpop.count_het()
            # ref = gt_subpop.count_hom_alt()
            alt = gt_subpop.count_hom_ref()
            print("{} hets, homalt: {} {}".format(pop, het, alt))
            chrstatdict[c][pop] = (seg, sing, doub)
        privalleles(ac_subpops, popdict)
        diff = fixdiff(ac_subpops, popdict)
    return(chrstatdict, diff)
