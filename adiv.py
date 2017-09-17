#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 18:15:38 2017

@author: scott
"""

import allel
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from jackknife import jackknife
sns.set_style('white')
sns.set_style('ticks')


def div_fx(ac_subpopscat, ac_subpopsdict, subpops, chrpos):
    """
    """
    poplist = ["Haiti", "Mali", "Kenya", "PNG"]
    for pop in poplist:
        pop_idx = subpops[pop]
        dip = (len(pop_idx)) * 2.0
        # prc missing
        prct = .80
        acu = allel.AlleleCountsArray(ac_subpopscat[pop][:])
        miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
        # seg = acu.is_segregating()
        # flt = miss * seg
        flt = miss
        blenW = int(round(sum(flt)/10))  # step size in SNP variants
        ac = allel.AlleleCountsArray(ac_subpopscat[pop].compress(flt,
                                     axis=0)[:, :2])
        # tajd
        d_window = allel.stats.diversity.moving_tajima_d(ac, blenW)
        d_blocks = allel.stats.moving_statistic(d_window, statistic=np.mean,
                                                size=1)
        d_m, d_se, d_z, d_val = jackknife(d_blocks)
        print("{}, TajD: {:.4f} +/- {:.4f}, z : {}".format(pop, d_m, d_se,
                                                           d_z))
        # pi
        pi = allel.stats.diversity.mean_pairwise_difference(ac)
        pi_blocks = allel.stats.moving_statistic(pi, statistic=np.mean,
                                                 size=blenW)
        pi_m, pi_se, pi_z, pi_val = jackknife(pi_blocks)
        print("{}, pi: {:.4f} +/- {:.4f}, z : {}".format(pop, pi_m, pi_se,
                                                         pi_z))
        theta = []
        for nchr in ac_subpopsdict.keys():
            acu = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop][:])
            miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
            # seg = acu.is_segregating()
            # flt = miss * seg
            flt = miss
            blenWc = int(round(sum(flt)/10))  # step size of variants
            acc = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop].
                                          compress(flt, axis=0)[:, :2])
            pos = chrpos[nchr][flt]
            # tajd
            d_window = allel.stats.diversity.moving_tajima_d(acc, blenWc)
            d_blocks = allel.stats.moving_statistic(d_window,
                                                    statistic=np.mean, size=1)
            d_m, d_se, d_z, d_val = jackknife(d_blocks)
            print("{} {}, TajD: {:.4f} +/- {:.4f}, z : {}".format(nchr, pop,
                                                                  d_m, d_se,
                                                                  d_z))
            # pi, already used for plotting similar to fst_hudson example
            pi = allel.stats.diversity.mean_pairwise_difference(acc)
            pi_blocks = allel.stats.moving_statistic(pi, statistic=np.mean,
                                                     size=blenWc)
            pi_m, pi_se, pi_z, pi_val = jackknife(pi_blocks)
            print("{} {}, pi: {:.4f} +/- {:.4f}, z : {}".format(nchr, pop,
                                                                pi_m, pi_se,
                                                                pi_z))
            # theta
            t = allel.stats.diversity.watterson_theta(pos, acc)
            theta.append(t)
            print("{} {}, theta: {:.6f}".format(nchr, pop, t))
            # theta windowed for plotting
            stop = max(pos)
            t_window = allel.stats.diversity.windowed_watterson_theta(
                    pos, acc, size=int(stop/1000), start=1, stop=stop)
            t_blocks = allel.stats.moving_statistic(t_window[0],
                                                    statistic=np.mean, size=10)
            t_m, t_se, t_z, t_val = jackknife(t_blocks)
            print("{} {}, theta: {:.6f} +/- {:.6f}, z : {}".format(nchr, pop,
                                                                   t_m, t_se,
                                                                   t_z))
        print("{}, theta: {:.6f} {:.6f}".format(pop, np.mean(theta),
                                                np.std(theta)))
    return(acc, pos)

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

