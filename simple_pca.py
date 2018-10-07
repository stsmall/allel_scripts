#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 16:24:36 2018

@author: scott
"""

from __future__ import division
from __future__ import print_function
import pandas as pd

# functions
from allel_class import Chr
import apca as apca
import autil as autil


chromlist = ['2L', '2R', '3L', '3R', 'X']
pops = ["Fun", "Like", "Long", "Par", "Van"]
pops = ["Ken", "Moz", "Uga", "Tan", "Zam", "Gha", "Kwa"]
for chrom in chromlist:
    var = Chr(pops, '../AnfunSG.liftover.biSNP.50.NoMiss.mac1.derived.noLinkSel.h5')
    meta = "../AnopSG.fun.info"
    meta = pd.read_csv(meta, delimiter=",")
    popdict = autil.subpops(var, meta, bypop=True, bykary=False)
    pop2color = autil.popcols(popdict)
    var.geno(chrom, meta)
    var.miss(var.gt, var.pos, 0)
    var.mac(var.gt, var.pos, 1)
    var.seg(var.gt, var.pos)
    # gn, thinp = autil.ldthin(var.gt, var.pos, "random")
    gn = var.gt.to_n_alt()[:]
    coords, model = apca.pca_fx(gn, meta, chrom, pop2color, False, pops)
