#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:40:01 2017

@author: scott
"""

import allel
import numpy as np


def jackknife(x):
    """
    """

    vals = np.empty(x.shape, dtype=float)
    x = np.ma.asarray(x)
    x.mask = np.zeros(x.shape, dtype=bool)
    for i in range(x.size):
        x.mask[i] = True
        vals[i] = np.mean(x)
        x.mask[i] = False
    n = x.size
    try:
        sv = ((n - 1) / n) * np.sum((vals - vals.mean()) ** 2)
    except ZeroDivisionError:
        se = 0.0000000001
    se = np.sqrt(sv)
    m = np.mean(vals)
    z = m / se
    # here the forumula is actually m - 0 / se, 0 is expected value
    return(m, se, z, vals)


def ldthin(geno, positions, method, mac=2, size=100, step=20, thresh=.1,
           iters=5):
    """.take if coord, .compress if mask
    """
    ac = geno.count_alleles()
    # only biallelic and mac or 2
    pca_selection = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > mac)

    print("Available site for PCA: {}".format(np.count_nonzero(pca_selection)))
    if 'random' in method:
        indices = np.nonzero(pca_selection)[0]
        indices_ds = np.random.choice(indices, size=50000, replace=False)
        indices_ds.sort()
        genotypes_pca = geno.take(indices_ds, axis=0)
        # sites with missing data can return error
        gn = genotypes_pca.to_n_alt()[:]
        pos = indices_ds
        print("{} Random SNPs selected for PCA".format(gn.shape[0]))
    else:
        genotypes_pca = geno.compress(pca_selection, axis=0)
        gn = genotypes_pca.to_n_alt()
        pos = positions[pca_selection]
        for i in range(iters):
            loc_unlinked = allel.locate_unlinked(gn, size=size, step=step,
                                                 threshold=thresh)
            n = np.count_nonzero(loc_unlinked)
            n_remove = gn.shape[0] - n
            print("iteration {} retaining {} removing {} variants".format(i+1,
                  n, n_remove))
            gn = gn.compress(loc_unlinked, axis=0)
            pos = positions[loc_unlinked]
    return(gn, allel.SortedIndex(pos))
