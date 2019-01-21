#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 16:34:02 2019

@author: scott
"""
from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib as mpl
# functions
from allel_class import Chr
import autil as autil
import adiv as av
import ald as ald
import adxy as adxy
import aplot as aplot
mpl.rcParams['pdf.fonttype'] = 42

meta = "AnopSG.55.info"
meta = pd.read_csv(meta, delimiter=",")
var = Chr('All', '2L.FSG.SNP.recode.h5')
popdict = autil.subpops(var, meta, bypop=True, bykary=False)
pop2color = autil.popcols(popdict)

# chrlist = np.unique(var.chrm[:])
chrlen = {}
with open("chr_info", 'r') as c:
    for line in c:
        x = line.strip().split()
        chrlen[x[0]] = int(x[1])

# RND
dxydict = {}
for c in chrlen.keys():
    var = Chr('All', '{}.FSG.SNP.recode.h5'.format(c))
    var.geno(c, meta)
    # var.miss(var.gt, var.pos, .20)
    # var.mac(var.gt, var.pos, 1)
    print("\nStats for Chromosome {}\n".format(c))
    # allele count object
    ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
    df_dxy = adxy.pairDxy(c, chrlen[c], ac_subpops, var.pos, plot=True)
    dxydict[c] = df_dxy

# Diversity statistics
pidict = {}
tajddict = {}
thetadict = {}
for c in chrlen.keys():
    var = Chr('All', '{}.FSG.SNP.recode.h5'.format(c))
    var.geno(c, meta)
    print("\nStats for Chromosome {}\n".format(c))
    # var.miss(var.gt, var.pos, .20)
    # var.mac(var.gt, var.pos, 1)
    # allele count object
    ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=1)
    pi = av.pi(c, chrlen[c], ac_subpops, var.pos, plot=True)
    pidict[c] = pi
    d = av.tajd(c, chrlen[c], ac_subpops, var.pos, plot=True)
    tajddict[c] = d
    t = av.theta(c, chrlen[c], ac_subpops, var.pos, plot=True)
    thetadict[c] = t

# Diversity boxplot; sumdict() autil then boxplot in aplot

theta = autil.catdict(thetadict)
aplot.divboxplot(theta, pop2color)
pi = autil.catdict(pidict)
aplot.divboxplot(pi, pop2color)
tajd = autil.catdict(tajddict)
aplot.divboxplot(tajd, pop2color)

# tajd histogram how do I get colors and transparency??
for pop in popdict.keys():
    b, bins, patches = mpl.pyplot.hist(tajd[pop], 50, density=True)

# LD decay plot
lddict = {}
for c in chrlen.keys():
    var = Chr('All', '{}.FSG.SNP.recode.h5'.format(c))
    var.geno(c, meta)
    print("\nStats for Chromosome {}\n".format(c))
    var.miss(var.gt, var.pos, .20)
    # allele count object
    ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=1)
    lddict[c] = ald.ld_decay(c, chrlen[c], ac_subpops, popdict,
                             pop2color, var)