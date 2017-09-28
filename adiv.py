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
from autil import jackknife
sns.set_style('white')
sns.set_style('ticks')


def div_plot(divdict, pops, chrom, chrsize, name, save=False):
    """
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    title = name
    for p in pops:
        nx = divdict[p][2][1]
        x = [(np.sum(i)-1)/2 for i in nx]  # need midpoints
        y = divdict[p][2][0]
        ax.plot(x, y, lw=.5, label=p)
    ax.set_ylabel(name)
    ax.set_xlabel('Chromosome {} position (bp)'.format(chrom))
    ax.set_xlim(0, chrsize)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    if save:
        fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


def pi(c, chrsize, ac_subpops, pos, plot=False, blenw=10000, nwindow=100):
    """
    """
    pidict = {}
    windlen = int(chrsize / nwindow)
    for x in ac_subpops.keys():
        acu = ac_subpops[x]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = pos[flt]
        ac = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                            axis=0)[:, :2])
        # pi
        pi = allel.windowed_diversity(posflt, ac, size=blenw)
        pi_m, pi_se, *f = jackknife(pi[0])
        pi_windowed = allel.windowed_diversity(posflt, ac, size=windlen,
                                               start=1, stop=chrsize)
        pidict[x] = (pi_m, pi_se, (pi_windowed[0], pi_windowed[1]))

    if plot:
        div_plot(pidict, list(ac_subpops.keys()), c, chrsize, "pi")
    return(pidict)


def theta(c, chrsize, ac_subpops, pos, plot=False, blenw=10000, nwindow=100):
    """
    """
    thetadict = {}
    windlen = int(chrsize / nwindow)
    for x in ac_subpops.keys():
        acu = ac_subpops[x]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = pos[flt]
        ac = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                            axis=0)[:, :2])
        # theta
        theta = allel.windowed_watterson_theta(posflt, ac, size=blenw)
        t_m, t_se, *t = jackknife(theta[0])
        theta_windowed = allel.windowed_watterson_theta(posflt, ac,
                                                        size=windlen, start=1,
                                                        stop=chrsize)
        thetadict[x] = (t_m, t_se, (theta_windowed[0], theta_windowed[1]))
    if plot:
        div_plot(thetadict, list(ac_subpops.keys()), c, chrsize, "theta")
    return(thetadict)


def tajd(c, chrsize, ac_subpops, pos, plot=False, blenw=10000, nwindow=100):
    """
    """
    tajddict = {}
    windlen = int(chrsize / nwindow)
    for x in ac_subpops.keys():
        acu = ac_subpops[x]
        flt = acu.is_segregating() & (acu.max_allele() == 1)
        print('retaining', np.count_nonzero(flt), 'SNPs')
        posflt = pos[flt]
        ac = allel.AlleleCountsArray(ac_subpops[x].compress(flt,
                                                            axis=0)[:, :2])
        # tajd
        tajd = allel.windowed_tajima_d(posflt, ac, size=blenw)
        d_m, d_se, *d = jackknife(tajd[0])
        tajd_windowed = allel.windowed_tajima_d(posflt, ac, size=windlen,
                                                start=1, stop=chrsize)
        # moving window of variants rather than based
#        tajd_sizevars = allel.moving_tajima_d(ac, size=size)
        tajddict[x] = (d_m, d_se, (tajd_windowed[0], tajd_windowed[1]))
    if plot:
        div_plot(tajddict, list(ac_subpops.keys()), c, chrsize, "Tajima's D")
    return(tajddict)
