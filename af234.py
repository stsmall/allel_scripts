#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:56:47 2017

@author: scott
"""

import allel
import numpy as np
# np.seterr(divide='ignore', invalid='ignore')
from itertools import combinations
from itertools import permutations
from autil import jackknife


def pF2(ac_subpops, blenw=10000):
    """
    """
    f2dict = {}
    # Patterson F2
    for x, y in combinations(ac_subpops.keys(), 2):
        acu = ac_subpops[x] + ac_subpops[y]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
#        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                             axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt,
                                                             axis=0)[:, :2])
        # Patterson f2
        f2 = allel.stats.admixture.patterson_f2(ac1, ac2)
        f2_blocks = allel.stats.moving_statistic(f2, statistic=np.nanmean,
                                                 size=blenw)
        f2_m, f2_se, fz, *f = jackknife(f2_blocks)
        f2dict["{}-{}".format(x, y)] = (f2_m, f2_se, fz)
    return(f2dict)


def pF3(ac_subpops, blenw=10000):
    """Patterson F3
    """
    f3dict = {}
    for x, y, z in permutations(ac_subpops.keys(), 3):
        acu = ac_subpops[x] + ac_subpops[y] + ac_subpops[z]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
#        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                             axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[y].compress(flt,
                                                             axis=0)[:, :2])
        ac3 = allel.AlleleCountsArray(ac_subpops[z].compress(flt,
                                                             axis=0)[:, :2])
#        f3, f3_se, f3z, *f = allel.stats.blockwise_patterson_f3(ac1, ac2, ac3,
#                                                                blenw)
        f3, f3_se, f3z, *f = allel.average_patterson_f3(ac1, ac2, ac3, blenw)
        f3dict["{}-{}-{}".format(x, y, z)] = (f3, f3_se, f3z)
    return(f3dict)


def pF4(ac_subpops, blenw=10000):
    """Patterson F4
    """
    f4dict = {}
    for w, x, y, z in permutations(ac_subpops.keys(), 4):
        acu = ac_subpops[w] + ac_subpops[x] + ac_subpops[y] + ac_subpops[z]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
#        posflt = pos[flt]
        ac1 = allel.AlleleCountsArray(ac_subpops[w].compress(flt,
                                                             axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                             axis=0)[:, :2])
        ac3 = allel.AlleleCountsArray(ac_subpops[y].compress(flt,
                                                             axis=0)[:, :2])
        ac4 = allel.AlleleCountsArray(ac_subpops[z].compress(flt,
                                                             axis=0)[:, :2])
        f4, f4_se, f4zscore, *f = allel.stats.blockwise_patterson_d(ac1, ac2,
                                                                    ac3, ac4,
                                                                    blenw)
        f4dict["{}-{}-{}-{}".format(w, x, y, z)] = (f4, f4_se, f4zscore)

    return(f4dict)
