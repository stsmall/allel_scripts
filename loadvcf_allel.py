#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 11:18:13 2017

@author: scott
"""

import h5py
import allel
import pandas as pd
from allel.py import popgen


def loaddata(h5, meta):
    """load h5 from file
    """
    anfun = popgen(anfun)
    # Prep and load callset
    callset_fn = h5
    anfun.call = h5py.File(callset_fn, mode='r')  # load
    anfun.var =
    anfun.gt = allel.GenotypeChunkedArray(anfun.call['calldata']['genotype'])
    anfun.ac = anfun.gt.count_alleles()[:]  # allele count
    flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)  # rmv singleton
    gf = g.compress(flt, axis=0)  # filtering
    gn = gf.to_n_alt()  # 2d matrix for correlation
    df_samples = pd.read_csv(meta, delimiter=' ')
    return(gn, df_samples)
