#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 13:47:33 2018

@author: scott
"""

from __future__ import division
from __future__ import print_function
import pandas as pd

# functions
from allel_class import Chr
import asfs as asfs
import autil as autil


chromlist = ['2L', '2R', '3L', '3R', 'X']
pops = ["Fun", "Like", "Long", "Par", "Van"]
for chrom in chromlist:
    var = Chr(pops, '/home/scott/Desktop/AnopSG_liftvcf/{}.SNP.recode.h5'.format(chrom))
    meta = "/home/scott/Desktop/AnopSG_liftvcf/AnopSG.55.info"
    meta = pd.read_csv(meta, delimiter=",")
    popdict = autil.subpops(var, meta, bypop=True, bykary=False)
    pop2color = autil.popcols(popdict)
    var.geno(chrom, meta)
    var.miss(var.gt, var.pos, 0)
    ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=1)
    sfs = asfs.sfs_plot(chrom, ac_subpops)
