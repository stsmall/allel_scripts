#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 16:55:12 2017

@author: scott
"""
import allel
import hmmlearn
from accessmasknp import access_mask_fx


def roh(subpops, chrpos, chrdict, nchr, chrmaskdict):
    """Runs of Homozygosity for each individual sample
    phet_roh: prob of obs a het within a RoH ~mutation rate
    phet_nonroh: >1 prob of obs het outside of RoH, ~nucleotide div
    """
    m = access_mask_fx(chrmaskdict, nchr)
    for pop in subpops:  # Haiti
        for indx in pop:  # 1, 2, 3, 4, ... 8
            gn = chrdict[nchr].genotypes[:, indx]
            het_mask = gn.is_het()
            gt = gn.compress(het_mask)
            call_mask = gn.is_called()  # all sites called
            posmask = call_mask * het_mask
            pos = chrdict[nchr].positions[posmask]
            df, prop = allel.stats.roh.roh_mhmm(gt, pos, phet_roh=0.001,
                                                phet_nonroh=(0.0025, 0.01),
                                                is_accessible=m)
    return(None)
